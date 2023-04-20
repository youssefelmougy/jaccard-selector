/******************************************************************
*
* An Asynchronous Distributed Actor-based Approach to Jaccard Similarity for Genome Comparisons
* 
 *****************************************************************/ 
/*! \file jaccard_kmer_selector.cpp
 * \brief Demo application that calculates jaccard similarity for k-mer sequence graph
 *          
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

#define THREADS shmem_n_pes()
#define MYTHREAD shmem_my_pe()

enum MailBoxType {RESPONSE, REQUEST};

typedef struct JaccardPkt {
    int64_t x;
    int64_t pos_row;
    int64_t pos_col;
    int64_t index_u;
} JaccardPkt;

typedef struct DegreePkt {
    int64_t i;
    int64_t src_idx;
} DegreePkt;

class JaccardSelector: public hclib::Selector<2, JaccardPkt> {
public:
    JaccardSelector(sparsemat_t* mat, sparsemat_t* intersection_mat) : mat_(mat), intersection_mat_(intersection_mat) {
        mb[REQUEST].process = [this] (JaccardPkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (JaccardPkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    sparsemat_t* mat_;
    sparsemat_t* intersection_mat_;

    void req_process(JaccardPkt pkg, int sender_rank) {
        JaccardPkt pkg2;
        for (int64_t uk = mat_->loffset[pkg.index_u]; uk < mat_->loffset[pkg.index_u+1]; uk++) {
            if (pkg.x == mat_->lnonzero[uk]) {
                pkg2.pos_row = pkg.pos_row;
                pkg2.pos_col = pkg.pos_col;
                send(RESPONSE, pkg2, pkg.pos_row % THREADS);
            }
        }
    }

    void resp_process(JaccardPkt pkg, int sender_rank) {
        int pos = 0;
        pkg.pos_row = pkg.pos_row / THREADS;
        for (int i = intersection_mat_->loffset[pkg.pos_row]; i < intersection_mat_->loffset[pkg.pos_row + 1]; i++) {
            if (pos == pkg.pos_col) {
                intersection_mat_->lvalue[i]++;
            }
            pos++;
        }
    }

};

class DegreeSelector: public hclib::Selector<2, DegreePkt> {
public:
    DegreeSelector(sparsemat_t* mat, sparsemat_t* kmer_mat) : mat_(mat), kmer_mat_(kmer_mat) {
        mb[REQUEST].process = [this] (DegreePkt pkt, int sender_rank) { 
            this->req_process(pkt, sender_rank);
        };
        mb[RESPONSE].process = [this] (DegreePkt pkt, int sender_rank) { 
            this->resp_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    sparsemat_t* mat_;
    sparsemat_t* kmer_mat_;

    void req_process(DegreePkt pkg, int sender_rank) {
        DegreePkt pkg2;
        pkg2.src_idx = pkg.src_idx;
        pkg2.i = kmer_mat_->loffset[pkg.i + 1] - kmer_mat_->loffset[pkg.i];
        send(RESPONSE, pkg2, sender_rank);
    }

    void resp_process(DegreePkt pkg, int sender_rank) {
        mat_->lnonzero[pkg.src_idx] = pkg.i;
    }

};

double get_edge_degrees(sparsemat_t* L, sparsemat_t* A2) {
    DegreeSelector* degreeSelector = new DegreeSelector(L, A2);

    hclib::finish([=]() {
        degreeSelector->start();
        DegreePkt degPKG;
    
        for (int64_t y = 0; y < L->lnumrows; y++) {
            for (int64_t k = L->loffset[y]; k < L->loffset[y+1]; k++) {
                degPKG.i = L->lnonzero[k] / THREADS;
                degPKG.src_idx = k;
                degreeSelector->send(REQUEST, degPKG, L->lnonzero[k] % THREADS);
            }
        }

        degreeSelector->done(REQUEST);
    });

    return 0;
}

double jaccard_selector(sparsemat_t* Jaccard_mat, sparsemat_t* A2) {
    // Start timing
    double t1 = wall_seconds();
    lgp_barrier();

    JaccardSelector* jacSelector = new JaccardSelector(A2, Jaccard_mat);

    get_edge_degrees(Jaccard_mat, A2);
    lgp_barrier();

    hclib::finish([=]() {
        jacSelector->start();

        JaccardPkt pkg;

        for (int64_t v = 0; v < A2->lnumrows; v++) {             //vertex v (local)
            for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
                int64_t v_nonzero = A2->lnonzero[k];                     //vertex u (possibly remote)
                int64_t row_num = v * THREADS + MYTHREAD;

                for (int64_t i_rows = row_num; i_rows < A2->numrows; i_rows++) {
                    // calculate intersection
                    pkg.index_u = i_rows / THREADS;
                    pkg.x = v_nonzero;
                    pkg.pos_row = i_rows;
                    pkg.pos_col = row_num;
                    jacSelector->send(REQUEST, pkg, i_rows % THREADS);
                }
            }
        }
        jacSelector->done(REQUEST);
        
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

