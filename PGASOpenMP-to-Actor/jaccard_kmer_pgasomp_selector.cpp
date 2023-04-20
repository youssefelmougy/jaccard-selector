/******************************************************************
*
* An Asynchronous Distributed Actor-based Approach to Jaccard Similarity for Genome Comparisons
* 
 *****************************************************************/ 
/*! \file jaccard_kmer_pgasomp_selector.cpp
 * \brief Demo application that calculates jaccard similarity for k-mer sequence graph
 *          THIS IS AN AUTOMATICALLY GENERATED ACTOR-BASED VERSION GIVEN A PGAS OPENMP PROGRAM (discussed in Section 6)
 *
 *      NOTE: this program is created for mxn k-mer matrices
 */

#include "/opt/cray/pe/perftools/23.02.0/include/pat_api.h"
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

    enum MailBoxType{__M1, __M2};
    Selector<2>* __selector = new Selector<2>(A2, Jaccard_mat);

    get_edge_degrees(Jaccard_mat, A2);
    lgp_barrier();

    hclib::finish([=]() {
        __selector->start();

        for (int64_t v = 0; v < A2->lnumrows; v++) {             //vertex v (local)
            for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
                int64_t v_nonzero = A2->lnonzero[k];                     //vertex u (possibly remote)
                int64_t row_num = v * shmem_n_pes() + shmem_my_pe();

                for (int64_t i_rows = row_num; i_rows < A2->numrows; i_rows++) {
                    // calculate intersection
                    int64_t __P1 = get_remote_pe(i_rows);
                    __selector->send(__M1, __P1, [=] () {
                        for (int64_t uk = A2->loffset[get_local_index(i_rows)]; uk < A2->loffset[get_local_index(i_rows + shmem_n_pes())]; uk++) {
                            if (v_nonzero == A2->lnonzero[uk]) {
                                int64_t __P2 = get_remote_pe(i_rows);
                                __selector->send(__M2, __P2, [=] () {
                                    int pos = 0;
                                    for (int i = Jaccard_mat->loffset[get_local_index(i_rows)]; i < Jaccard_mat->loffset[get_local_index(i_rows + shmem_n_pes())]; i++) {
                                        if (pos == row_num) {
                                            Jaccard_mat->lvalue[i]++;
                                        }
                                        pos++;
                                    }
                                });
                            }
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

        T0_fprintf(stderr,"Running jaccard on %d threads\n", THREADS);

        // Read/create k-mer matrix
        sparsemat_t* A2 = read_matrix_mm_to_dist("kmer_matrix.mtx");
        A2 = transpose_matrix(A2);
        lgp_barrier();


        // Initialize Jaccard similarity matrix
        sparsemat_t* Jaccard_mat = random_graph(5000, FLAT, DIRECTED_WEIGHTED, LOOPS, 1, 12345); if (MYTHREAD==1) Jaccard_mat->lnonzero[0] = 0; tril(Jaccard_mat, -1); for (int r = 0; r < Jaccard_mat->lnumrows; r++) {for (int rr = Jaccard_mat->loffset[r]; rr < Jaccard_mat->loffset[r+1]; rr++) {Jaccard_mat->lvalue[rr] = 0.0;}}
        lgp_barrier();

        T0_fprintf(stderr, "K-mer Matrix is %ldx%ld and has %ld nonzeros.\n", A2->numcols, A2->numrows, A2->nnz);
        T0_fprintf(stderr, "\nJaccard Similarity Matrix is %ldx%ld and has %ld values.\n", Jaccard_mat->numrows, Jaccard_mat->numrows, Jaccard_mat->nnz);

        T0_fprintf(stderr, "\nRunning Jaccard Similarity K-mers (selector): \n");

        double laptime_jaccard = 0.0;
        
        PAT_region_begin(1, "selector_jaccard"); // region begin for HWPC evaluation using CrayPat

        // Running selector model for jaccard similarity on k-mer matrix
        laptime_jaccard = jaccard_selector(Jaccard_mat, A2);
        lgp_barrier();

        PAT_region_end(1); // region end for HWPC evaluation using CrayPat

        T0_fprintf(stderr, " %8.5lf seconds\n", laptime_jaccard);
        lgp_barrier();

        lgp_finalize();
        
    });

    return 0;
}

