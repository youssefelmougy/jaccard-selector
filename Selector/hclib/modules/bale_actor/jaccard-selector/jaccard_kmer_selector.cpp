/******************************************************************
*
* Asynchronous Distributed Actor-based Approach to Jaccard Similarity for Genome Comparisons
* 
 *****************************************************************/ 
/*! \file jaccard_kmer_selector.cpp
 * \brief Demo application that calculates jaccard similarity for k-mer genome sequences 
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

int64_t binary_search(int64_t l_, int64_t r_, int64_t val_, sparsemat_t* mat_, int64_t flag_) { // custom binary search function
    while (l_ <= r_) {
        int m_ = l_ + (r_ - l_) / 2;
        int64_t comp_var_;
        if (flag_) { comp_var_ = m_ - l_; } else { comp_var_ = mat_->lnonzero[m_]; }
        if (comp_var_ == val_) { return m_; }
        if (comp_var_ < val_) { l_ = m_ + 1; } else { r_ = m_ - 1; }
    }
    return -1;
}

class JaccardSelector: public hclib::Selector<1, JaccardPkt> {
public:
    JaccardSelector(sparsemat_t* mat, sparsemat_t* intersection_mat) : mat_(mat), intersection_mat_(intersection_mat) {
	    mb[0].process = [this] (JaccardPkt pkt, int sender_rank) { 
	        this->req_process(pkt, sender_rank);
        };
    }

private:
    //shared variables
    sparsemat_t* mat_;
    sparsemat_t* intersection_mat_;

    void req_process(JaccardPkt pkg, int sender_rank) {
        JaccardPkt pkg2;
        int64_t search_ = binary_search(mat_->loffset[pkg.index_u], mat_->loffset[pkg.index_u + 1] - 1, pkg.x, mat_, 0);
        if (search_ != -1) {
            pkg.pos_row = pkg.pos_row / THREADS;
            int64_t search__ = binary_search(intersection_mat_->loffset[pkg.pos_row], intersection_mat_->loffset[pkg.pos_row + 1] - 1, pkg.pos_col, mat_, 1);
            if (search__ != -1) { intersection_mat_->lvalue[search__]++; }
        }
    }
};

class DegreeSelector: public hclib::Selector<2, DegreePkt> {
public:
    DegreeSelector(sparsemat_t* mat, sparsemat_t* kmer_mat) : mat_(mat), kmer_mat_(kmer_mat) {
	    mb[0].process = [this] (DegreePkt pkt, int sender_rank) { 
	        this->req_process(pkt, sender_rank);
        };
        mb[1].process = [this] (DegreePkt pkt, int sender_rank) { 
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
        send(1, pkg2, sender_rank); 
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
                degreeSelector->send(0, degPKG, L->lnonzero[k] % THREADS); 
            }
        }

        degreeSelector->done(0); 
    });

    return 0;
}

double jaccard_selector(sparsemat_t* Jaccard_mat, sparsemat_t* A2) {
    double t1 = wall_seconds();
    lgp_barrier();

    JaccardSelector* jacSelector = new JaccardSelector(A2, Jaccard_mat);

    get_edge_degrees(Jaccard_mat, A2);
    lgp_barrier();

    hclib::finish([=]() {
        jacSelector->start();

        JaccardPkt pkg;

        for (int64_t v = 0; v < A2->lnumrows; v++) {         	         //vertex v (local)
            for (int64_t k = A2->loffset[v]; k < A2->loffset[v+1]; k++) {
                int64_t v_nonzero = A2->lnonzero[k];                     //vertex u (possibly remote)
                int64_t row_num = v * THREADS + MYTHREAD;

                for (int64_t i_rows = row_num; i_rows < A2->numrows; i_rows++) {
                    // calculate intersection
                    pkg.index_u = i_rows / THREADS;
                    pkg.x = v_nonzero;
                    pkg.pos_row = i_rows;
                    pkg.pos_col = row_num;
                    jacSelector->send(0, pkg, i_rows % THREADS);
                }
            }
        }
        jacSelector->done(0);
    });

    lgp_barrier();

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
        
        PAT_region_begin(1, "selector_jaccard"); // region begin for HWPC evaluation using CrayPat

        // Running selector model for jaccard similarity on k-mer matrix
        laptime_jaccard = jaccard_selector(Jaccard_mat, A2);
        lgp_barrier();

        PAT_region_end(1); // region begin for HWPC evaluation using CrayPat

        //T0_fprintf(stderr, "\nJACCARD KMER MATRIX: \n"); lgp_barrier();
        //print_matrix(A2);
        //printf("\n");

        //T0_fprintf(stderr, "\nJACCARD SIMILARITY MATRIX: \n"); lgp_barrier();
        //print_matrix(Jaccard_mat);
        //printf("\n");    

        T0_fprintf(stderr, " %8.5lf seconds\n", laptime_jaccard);
        lgp_barrier();

        lgp_finalize();
    });

    return 0;
}
