/******************************************************************
*
* Asynchronous Distributed Actor-based Approach to Jaccard Similarity for Genome Comparisons
* 
 *****************************************************************/ 
/*! \file jaccard_kmer_pgasomp_selector.cpp
 * \brief Demo application that calculates jaccard similarity for k-mer sequence graph
 *          THIS IS AN AUTOMATICALLY GENERATED ACTOR-BASED VERSION GIVEN A PGAS OPENMP PROGRAM (discussed in Section 5)
 *
 *      NOTE: this program is created for mxn k-mer matrices
 */

#include <math.h>
#include <shmem.h>
extern "C" {
#include "spmat.h"
}
#include "selector.h"

double get_edge_degrees(sparsemat_t* L, sparsemat_t* A2) {
    enum MailBoxType{__M1, __M2};
    Selector<2>* __selector = new Selector<2>(L, A2);

    hclib::finish([=]() {
        __selector->start();

        int64_t __P1 = shmem_my_pe();
        for (int64_t y = 0; y < L->lnumrows; y++) {
            for (int64_t k = L->loffset[y]; k < L->loffset[y+1]; k++) {
                int64_t __P2 = get_remote_pe(L->lnonzero[k]);
                __selector->send(__M1, __P2, [=] () {
                    int64_t __t1 = A2->loffset[get_local_index(L->lnonzero[k] + shmem_n_pes())] - A2->loffset[get_local_index(L->lnonzero[k])];
                    __selector->send(__M2, __P1, [=] () {
                        L->lnonzero[k] = __t1;
                    });
                });
            }
        }
        __selector->done(__M1);
    });

    return 0;
}

double jaccard_selector(sparsemat_t* Jaccard_mat, sparsemat_t* A2) {
    // Start timing
    double t1 = wall_seconds();
    lgp_barrier();

    enum MailBoxType{__M1};
    Selector<1>* __selector = new Selector<1>(A2, Jaccard_mat);

    get_edge_degrees(Jaccard_mat, A2);
    lgp_barrier();

    hclib::finish([=]() {
        __selector->start();

        for (int64_t v = 0; v < A2->lnumrows; v++) {             //vertex v (local)
            int64_t v_global = v * shmem_n_pes() + shmem_my_pe();
            for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
                int64_t u = A2->lnonzero[k];                     //vertex u (possibly remote)
                for (int64_t i_rows = v_global; i_rows < A2->numrows; i_rows++) {
                   // calculate intersection
                    int64_t __P1 = get_remote_pe(i_rows);
                    __selector->send(__M1, __P1, [u, v_global, i_rows] () {
                        if (binary_search(A->lnonzero[A->loffset[get_local_index(i_rows)]:A->loffset[get_local_index(i_rows+shmem_n_pes())]], u)) {
                            int64_t pos = binary_search(Jaccard_mat->loffset[get_local_index(i_rows)], Jaccard_mat->loffset[get_local_index(i_rows + shmem_n_pes())], v_global);
                            Jaccard_mat->lvalue[pos]++;
                        }
                    });
                }
            }
        }
        __selector->done(__M1);   
    });

    lgp_barrier();

    // calculate Jaccard similarity
    for (int64_t v = 0; v < Jaccard_mat->lnumrows; v++) {
        int64_t deg_v = A2->loffset[v+1] - A2->loffset[v];
        for (int64_t k = Jaccard_mat->loffset[v]; k < Jaccard_mat->loffset[v+1]; k++) {
            Jaccard_mat->lvalue[k] = (double)Jaccard_mat->lvalue[k] / ((double)Jaccard_mat->lnonzero[k] + (double)deg_v - (double)Jaccard_mat->lvalue[k]);
        }
    }

    lgp_barrier();
    t1 = wall_seconds() - t1;
    return t1;
}


int main(int argc, char* argv[]) {
    const char *deps[] = { "system", "bale_actor" };
    hclib::launch(deps, 2, [=] {
		    
        double time_main_total = 0.0;
        time_main_total = wall_seconds();

        T0_fprintf(stderr,"Running jaccard on %d threads\n", THREADS);

        // Read in variables
        int opt;
        char filename[64];
        int64_t mtx_num_cols;
        int64_t mtx_num_rows;
        double nonzero_prob_kmer;
        int64_t read_graph = 0L;           // read graph from a file
        while ((opt = getopt(argc, argv, "m:n:p:f:")) != -1) {
            switch (opt) {
                case 'm': sscanf(optarg,"%ld", &mtx_num_rows);  break;
                case 'n': sscanf(optarg,"%ld", &mtx_num_cols); break;
                case 'p': sscanf(optarg,"%lg", &nonzero_prob_kmer); break;
                case 'f': read_graph = 1; sscanf(optarg,"%s", filename); break;
                default:  break;
            }
        }

        // Read/Generate jaccard k-mer matrix
        sparsemat_t* A2;
        if (read_graph) {
            A2 = read_matrix_mm_to_dist(filename);
        } else {
            A2 = generate_kmer_matrix(mtx_num_rows, mtx_num_cols, nonzero_prob_kmer);
        }
        A2 = transpose_matrix(A2);
        mtx_num_rows = A2->numrows;
        mtx_num_cols = A2->numcols;
        lgp_barrier();

        T0_fprintf(stderr, "K-mer Matrix is %ldx%ld and has %ld nonzeros.\n\n", A2->numcols, A2->numrows, A2->nnz);

        // Generate jaccard similarity matrix
        sparsemat_t * Jaccard_mat = generate_dense_tril_matrix(mtx_num_rows, 1); lgp_barrier();
        T0_fprintf(stderr, "Jaccard Similarity Matrix is %ldx%ld and has %ld values.\n", Jaccard_mat->numrows, Jaccard_mat->numrows, Jaccard_mat->nnz);

        T0_fprintf(stderr, "\nRunning Jaccard Similarity K-mers (selector): \n");

        double laptime_jaccard = 0.0;

        // Running selector model for jaccard similarity on k-mer matrix
        laptime_jaccard = jaccard_selector(Jaccard_mat, A2);
        lgp_barrier();

        T0_fprintf(stderr, " %8.5lf seconds\n", laptime_jaccard);
        lgp_barrier();

        lgp_finalize();
    });

    return 0;
}

