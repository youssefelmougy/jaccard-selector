/*
 * This file uses boost library
 * To compile: mpicxx -std=c++0x -O0 -g -Wall -Wno-format -DPMPI -DMPIIO -o fasta_reader fasta_reader.cxx -I/../ctf/include  -L/../ctf/lib -lctf -lblas -L/usr/local/opt/openblas/lib -L/../ctf/scalapack/build/lib -L/usr/local/gfortran/lib -llapack -lblas -lscalapack -lgfortran -lz -lboost_iostreams
 *
 */
#include <ctf.hpp>
#include <string>
#include <cstring>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


using namespace CTF;

void process_fasta_files(int64_t m, int64_t n, int64_t k, char *infolderPath, char *outfolderPath, const char *listfile, double *perc, int reverse_complement, World & dw)
{
  // nfiles: number of files this MPI process handles
  int64_t nfiles;
  nfiles = (n / dw.np) + (dw.rank < (n % dw.np));
  int64_t maxfiles;
  // max files are handled by rank 0
  // variable used to sync A.write()s across processes
  maxfiles = (n / dw.np) + (0 < (n % dw.np));
  uint64_t maxkmer = (uint64_t)1 << ((k - 1) * 2);
  maxkmer = *perc * maxkmer;
  // std::cout << maxkmer << endl;

  std::map<char, int> atgc = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}};
  // reverse_complement A->T, C->G, G->C, T->A
  std::map<char, int> atgc_reverse = {{'A', 3}, {'C', 2}, {'G', 1}, {'T', 0}};
  FILE *fplist;
  if (listfile != nullptr) {
    fplist = fopen(listfile, "r");
    if (fplist == nullptr && dw.rank == 0) {
      printf("I am unable to open file: %s\n", listfile);
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
    char dummy[9000];
    // Each rank starts from its corresponding file
    for (int64_t i = 0; i < dw.rank; i++) {
      fscanf(fplist, "%s", dummy);
    }
  }
  for (int64_t i = 0; i < maxfiles; i++) {
    if (i >= nfiles) {
      continue;
    }

    std::set<uint64_t> kmers;
    if (fplist != nullptr) {
      // Read files in lexicographic order
      char dummy[9000];
      fscanf(fplist, "%s", dummy);
      std::string ss = std::string(infolderPath) + "/" + std::string(dummy);
      std::ifstream file(ss.c_str(), std::ios_base::in | std::ios_base::binary);
      if (!file.is_open()) {
        printf("I am rank: %d, I was unable to open file: %s\n", dw.rank, ss.c_str());
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      // printf("I am rank: %d, I am opening file: %s\n", dw.rank, ss.c_str());
      boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
      inbuf.push(boost::iostreams::gzip_decompressor());
      inbuf.push(file);
      //Convert streambuf to istream
      std::istream instream(&inbuf);

      ss = std::string(outfolderPath) + "/" + string(dummy).substr(0, string(dummy).size() - 9) + ".txt";
      std::ofstream wfile(ss.c_str());
      if (!wfile.is_open()) {
        printf("I am rank: %d, I was unable to open file: %s\n", dw.rank, ss.c_str());
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, -1);
      }
      // printf("I am rank: %d, I am opening file: %s\n", dw.rank, ss.c_str());

      //Iterate lines
      std::string line;
      std::string pline;
      uint64_t ik = 0;
      uint64_t kmerv = 0;
      // can do away without a base; but the code is easier to understand this way
      uint64_t kmerv_b = 1;
      std::map<char, int>::iterator it;
      std::string::size_type il;
      bool process_rev_complement = false;
      while (std::getline(instream, line)) {
        // std::cout << line << std::endl;
        if (line.find(">") != std::string::npos) {
          // header found, process a new read
          ik = 0;
          kmerv = 0; kmerv_b = 1;
          if (reverse_complement) {
            // process reverse_complement; the string can be from two separate lines without the delimiter 
            if (pline.size() == 0) continue;
            // std::reverse(pline.begin(), pline.end());
            line = pline;
            std::cout << "reverse: " << line << std::endl;
            process_rev_complement = true;
          }
          else {
            continue;
          }
        }
        if (ik != 0) {
          // new line, but not separated by a header, so the read continues
          // prefix the previous line
          line = pline + line;
          il = pline.size();
        }
        else {
          il = 0;
        }
        for(; il < line.size(); il++) {
          if (!process_rev_complement) {
            it = atgc.find(toupper(line[il]));
          }
          else {
            it = atgc_reverse.find(toupper(line[il]));
          }
          if ((it != atgc.end() && !process_rev_complement) || (it !=atgc_reverse.end() && process_rev_complement)) {
            kmerv += kmerv_b * it->second;
            kmerv_b *= 4;
            //printf("kmerv: %lld char: %c val: %d\n", kmerv, line[il], it->second);
            ik++;
            if (ik == k) {
              // store the k-mer value
              // std::cout << kmerv << endl;
              if (kmerv <= maxkmer) {
                // wfile << kmerv << "\n";
                kmers.insert(kmerv);
              }
              ik = 0;
              kmerv = 0; kmerv_b = 1;
              il = il - (k - 1);
              assert (il >= 0);
            }
          }
          else {
            // found a non-atgc character
            ik = 0;
            kmerv = 0; kmerv_b = 1;
          }
        }
        // store this line as the previous line
        pline = line;
        if (process_rev_complement == true) {
          process_rev_complement = false;
          ik = 0;
          kmerv = 0; kmerv_b = 1;
        }
      }
      // printf("kmers.size(): %lld\n", kmers.size());
      // Maxkmer, and number of kmers are part of the file header
      wfile << maxkmer << " " << kmers.size() << "\n";
      // write the kmers
      for(auto const& value: kmers) {
        wfile << value << "\n";
      }
      kmers.clear();
      file.close();
      wfile.close();
      for (int64_t i = 0; i < (dw.np - 1); i++) {
        fscanf(fplist, "%s", dummy);
      }
    }
  }
} 

char* getCmdOption(char ** begin,
                   char ** end,
                   const   std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end){
    return *itr;
  }
  return 0;
}

int main(int argc, char ** argv){
  int rank, np;
  int64_t m, n, k;
  double perc;
  int const in_num = argc;
  char ** input_str = argv;
  char *infolderPath = nullptr;
  char *outfolderPath = nullptr;
  char *listfile = nullptr;
  int reverse_complement;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  {
    World dw(MPI_COMM_WORLD);
    
    if (getCmdOption(input_str, input_str+in_num, "-m")){
      m = atoll(getCmdOption(input_str, input_str+in_num, "-m"));
      if (m < 0) m = 1023;
    } else m = 1023;

    if (getCmdOption(input_str, input_str+in_num, "-n")){
      n = atoll(getCmdOption(input_str, input_str+in_num, "-n"));
      if (n < 0) n = 2;
    } else n = 2;
    
    if (getCmdOption(input_str, input_str+in_num, "-k")){
      k = atoll(getCmdOption(input_str, input_str+in_num, "-k"));
      if (k < 0) k = 3;
    } else k = 3;

    if (getCmdOption(input_str, input_str+in_num, "-perc")){
      perc = atof(getCmdOption(input_str, input_str+in_num, "-perc"));
      if (perc < 0) perc = .1;
    } else perc = .1;
 
    if (getCmdOption(input_str, input_str+in_num, "-lfile")) {
       listfile = getCmdOption(input_str, input_str+in_num, "-lfile");
     } else listfile = nullptr;
    
    if (getCmdOption(input_str, input_str+in_num, "-infolderPath")) {
       infolderPath = getCmdOption(input_str, input_str+in_num, "-infolderPath");
     } else infolderPath = nullptr;

    if (getCmdOption(input_str, input_str+in_num, "-outfolderPath")) {
       outfolderPath = getCmdOption(input_str, input_str+in_num, "-outfolderPath");
     } else outfolderPath = nullptr;
    
    if (getCmdOption(input_str, input_str+in_num, "-reverse_complement")){
      reverse_complement = atoi(getCmdOption(input_str, input_str+in_num, "-reverse_complement"));
      if (reverse_complement < 0) reverse_complement = 0;
    } else reverse_complement = 0;
    
    if ((listfile == nullptr || infolderPath == nullptr || outfolderPath == nullptr) && rank == 0) {
      printf("Error in the command line parameters");
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, -2);
    }
    double stime;
    double etime;
    stime = MPI_Wtime();
    process_fasta_files(m, n, k, infolderPath, outfolderPath, listfile, &perc, reverse_complement, dw);
    etime = MPI_Wtime() - stime;
    std::cout << "Total time taken: " << etime << endl;
  }
  MPI_Finalize();
}

