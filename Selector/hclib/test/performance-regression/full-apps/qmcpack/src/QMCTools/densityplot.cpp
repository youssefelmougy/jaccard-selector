#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include "OhmmsData/OhmmsElementBase.h"
#include "OhmmsData/HDFAttribIO.h"
#include "OhmmsPETE/TinyVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
using namespace qmcplusplus;
using namespace std;

int print_help()
{
  cout << "Usage: density [--fileroot string | --diff root1 root2 ] --np int " << endl;
  return 1;
}

struct DensityObserver
{
  typedef TinyVector<double,3> PosType;
  int npts;
  int nsamples;
  double rmin;
  double rmax;
  double delta;
  double deltainv;
  Matrix<double> value;
  DensityObserver():rmin(-7.92),rmax(7.92),delta(0.16),nsamples(0)
  {
    deltainv=1.0/delta;
    npts=static_cast<int>((rmax-rmin)/delta)+1;
    value.resize(npts,npts);
    value=0.0;
    cout << " Density grid " << npts << endl;
  }

  inline void accumulate(const vector<PosType>& pos, int nat)
  {
    nsamples+=pos.size();
    for(int i=0; i<nat; i++)
    {
      double x=pos[i][0];
      double y=pos[i][1];
      if(x>rmin && x<rmax && y >rmin && y<rmax)
      {
        value(static_cast<int>((x-rmin)*deltainv),static_cast<int>((y-rmin)*deltainv))+=1.0;
      }
    }
  }

  void getData(const char* fname, int nproc);

  void print(const char* fname)
  {
    ofstream fout(fname);
    fout.setf(ios::scientific, ios::floatfield);
    fout.setf(ios::left,ios::adjustfield);
    fout.precision(6);
    double norm=1.0/static_cast<double>(nsamples);
    double x=rmin;
    for(int ix=0; ix<npts; ix++,x+=delta)
    {
      double y=rmin;
      for(int jx=0; jx<npts; jx++,y+=delta)
        fout << setw(20) << y << setw(20) << x << setw(20) << norm*value(ix,jx) << endl;
      fout << endl;
    }
  }

};

int main(int argc, char** argv)
{
  if(argc<2)
    return print_help();
  int iargc=0;
  int nproc=1;
  vector<string> fnamein(2);
  string h5fileroot("0");
  bool getdiff=false;
  while(iargc<argc)
  {
    string c(argv[iargc]);
    if(c == "--fileroot")
    {
      h5fileroot=argv[++iargc];
    }
    else
      if(c == "--diff")
      {
        fnamein[0]=argv[++iargc];
        fnamein[1]=argv[++iargc];
        getdiff=true;
      }
      else
        if(c == "--np")
        {
          nproc=atoi(argv[++iargc]);
        }
        else
          if(c == "--help" || c == "-h")
          {
            return print_help();
          }
    ++iargc;
  }
  if(getdiff)
  {
    cout << "Going to compare two density " <<  fnamein[0] << " " << fnamein[1] << endl;
    DensityObserver recorder1;
    recorder1.getData(fnamein[0].c_str(),nproc);
    DensityObserver recorder2;
    recorder2.getData(fnamein[1].c_str(),nproc);
    recorder1.value -= recorder2.value;
    recorder1.print("rho_diff.dat");
  }
  else
  {
    DensityObserver recorder;
    recorder.getData(h5fileroot.c_str(),nproc);
  }
  return 0;
}

void DensityObserver::getData(const char* froot, int nproc)
{
  typedef TinyVector<double,3> PosType;
  hsize_t dimIn[3],dimTot[3];
  char fname[256],coordname[128];
  for(int ip=0; ip<nproc; ip++)
  {
    if(nproc>1)
      sprintf(fname,"%s.p%03d.config.h5",froot,ip);
    else
      sprintf(fname,"%s.config.h5",froot);
    hid_t fileID =  H5Fopen(fname,H5F_ACC_RDWR,H5P_DEFAULT);
    int nconf=1;
    hid_t mastercf = H5Gopen(fileID,"config_collection");
    hid_t h1=H5Dopen(mastercf,"NumOfConfigurations");
    H5Dread(h1, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&(nconf));
    H5Dclose(h1);
    cout << "Adding data in " << fname << " number of configurations = " << nconf << endl;
    for(int iconf=0; iconf<nconf; iconf++)
    {
      sprintf(coordname,"config%04d/coord",iconf);
      hid_t dataset = H5Dopen(mastercf,coordname);
      hid_t dataspace = H5Dget_space(dataset);
      int rank = H5Sget_simple_extent_ndims(dataspace);
      int status_n = H5Sget_simple_extent_dims(dataspace, dimTot, NULL);
      vector<PosType> pos(dimTot[0]*dimTot[1]);
      hid_t ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(pos[0][0]));
      H5Dclose(dataset);
      H5Sclose(dataspace);
      accumulate(pos,dimTot[0]*dimTot[1]);
    }
  }
  string outfname(froot);
  outfname.append(".den.dat");
  print(outfname.c_str());
}
