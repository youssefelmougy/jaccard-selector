#include <stdio.h>
__device__ inline int hclib_get_current_worker() {
    return blockIdx.x * blockDim.x + threadIdx.x;
}

template<class functor_type>
__global__ void wrapper_kernel(unsigned iter_offset, unsigned niters, functor_type functor) {
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < niters) {
        functor(iter_offset + tid);
    }
}
template<class functor_type>
static void kernel_launcher(const char *kernel_lbl, unsigned iter_offset, unsigned niters, functor_type functor) {
    const int threads_per_block = 256;
    const int nblocks = (niters + threads_per_block - 1) / threads_per_block;
    functor.transfer_to_device();
    const unsigned long long start = capp_current_time_ns();
    wrapper_kernel<<<nblocks, threads_per_block>>>(iter_offset, niters, functor);
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "CUDA Error while synchronizing kernel - %s\n", cudaGetErrorString(err));
        exit(2);
    }
    const unsigned long long end = capp_current_time_ns();
    fprintf(stderr, "%s %llu ns\n", kernel_lbl, end - start);
    functor.transfer_from_device();
}
#ifdef __cplusplus
#ifdef __CUDACC__
#endif
#endif
#include <stdio.h>
#include <omp.h>

extern int omp_num_threads;

#define BS 16

#define AA(_i,_j) a[offset*size+_i*size+_j+offset]
#define BB(_i,_j) a[_i*size+_j]

void lud_diagonal_omp (float* a, int size, int offset)
{
    int i, j, k;
    for (i = 0; i < BS; i++) {

        for (j = i; j < BS; j++) {
            for (k = 0; k < i ; k++) {
                AA(i,j) = AA(i,j) - AA(i,k) * AA(k,j);
            }
        }
   
        float temp = 1.f/AA(i,i);
        for (j = i+1; j < BS; j++) {
            for (k = 0; k < i ; k++) {
                AA(j,i) = AA(j,i) - AA(j,k) * AA(k,i);
            }
            AA(j,i) = AA(j,i)*temp;
        }
    }

}

// implements block LU factorization 
class pragma52_omp_parallel_hclib_async {
    private:
        void **host_allocations;
        size_t *host_allocation_sizes;
        unsigned nallocations;
        void **device_allocations;
    float* a;
    float* h_a;
    int size;
    volatile int offset;
    int chunk_idx;

    public:
        pragma52_omp_parallel_hclib_async(float* set_a,
                int set_size,
                int set_offset,
                int set_chunk_idx) {
            h_a = set_a;
            size = set_size;
            offset = set_offset;
            chunk_idx = set_chunk_idx;

        }

    void transfer_to_device() {
        int i;
        cudaError_t err;

        a = NULL;

        get_underlying_allocations(&host_allocations, &host_allocation_sizes, &nallocations, 1, h_a);
        device_allocations = (void **)malloc(nallocations * sizeof(void *));
        for (i = 0; i < nallocations; i++) {
            err = cudaMalloc((void **)&device_allocations[i], host_allocation_sizes[i]);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
            err = cudaMemcpy((void *)device_allocations[i], (void *)host_allocations[i], host_allocation_sizes[i], cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
            if (a == NULL && (char *)h_a >= (char *)host_allocations[i] && ((char *)h_a - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_a - (char *)host_allocations[i]);
                memcpy((void *)(&a), (void *)(&tmp), sizeof(void *));
            }
        }

        assert(a || h_a == NULL);

    }

    void transfer_from_device() {
        cudaError_t err;
        int i;
        for (i = 0; i < nallocations; i++) {
            err = cudaMemcpy((void *)host_allocations[i], (void *)device_allocations[i], host_allocation_sizes[i], cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
            err = cudaFree(device_allocations[i]);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
        }
    }

        __device__ void operator()(int chunk_idx) {
            for (int __dummy_iter = 0; __dummy_iter < 1; __dummy_iter++) {
                {
            int i, j, k, i_global, j_global, i_here, j_here;
            float sum;           
            float temp[BS*BS] __attribute__ ((aligned (64)));

            for (i = 0; i < BS; i++) {
for (j =0; j < BS; j++){
                    temp[i*BS + j] = a[size*(i + offset) + offset + j ];
                }
            }
            i_global = offset;
            j_global = offset;
            
            // processing top perimeter
            //
            j_global += BS * (chunk_idx+1);
            for (j = 0; j < BS; j++) {
                for (i = 0; i < BS; i++) {
                    sum = 0.f;
                    for (k=0; k < i; k++) {
                        sum += temp[BS*i +k] * BB((i_global+k),(j_global+j));
                    }
                    i_here = i_global + i;
                    j_here = j_global + j;
                    BB(i_here, j_here) = BB(i_here,j_here) - sum;
                }
            }

            // processing left perimeter
            //
            j_global = offset;
            i_global += BS * (chunk_idx + 1);
            for (i = 0; i < BS; i++) {
                for (j = 0; j < BS; j++) {
                    sum = 0.f;
                    for (k=0; k < j; k++) {
                        sum += BB((i_global+i),(j_global+k)) * temp[BS*k + j];
                    }
                    i_here = i_global + i;
                    j_here = j_global + j;
                    a[size*i_here + j_here] = ( a[size*i_here+j_here] - sum ) / a[size*(offset+j) + offset+j];
                }
            }

        }
            }
        }
};

class pragma106_omp_parallel_hclib_async {
    private:
        void **host_allocations;
        size_t *host_allocation_sizes;
        unsigned nallocations;
        void **device_allocations;
    volatile int offset;
    int chunk_idx;
    volatile int chunks_in_inter_row;
    float* a;
    float* h_a;
    int size;

    public:
        pragma106_omp_parallel_hclib_async(int set_offset,
                int set_chunk_idx,
                int set_chunks_in_inter_row,
                float* set_a,
                int set_size) {
            offset = set_offset;
            chunk_idx = set_chunk_idx;
            chunks_in_inter_row = set_chunks_in_inter_row;
            h_a = set_a;
            size = set_size;

        }

    void transfer_to_device() {
        int i;
        cudaError_t err;

        a = NULL;

        get_underlying_allocations(&host_allocations, &host_allocation_sizes, &nallocations, 1, h_a);
        device_allocations = (void **)malloc(nallocations * sizeof(void *));
        for (i = 0; i < nallocations; i++) {
            err = cudaMalloc((void **)&device_allocations[i], host_allocation_sizes[i]);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
            err = cudaMemcpy((void *)device_allocations[i], (void *)host_allocations[i], host_allocation_sizes[i], cudaMemcpyHostToDevice);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
            if (a == NULL && (char *)h_a >= (char *)host_allocations[i] && ((char *)h_a - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_a - (char *)host_allocations[i]);
                memcpy((void *)(&a), (void *)(&tmp), sizeof(void *));
            }
        }

        assert(a || h_a == NULL);

    }

    void transfer_from_device() {
        cudaError_t err;
        int i;
        for (i = 0; i < nallocations; i++) {
            err = cudaMemcpy((void *)host_allocations[i], (void *)device_allocations[i], host_allocation_sizes[i], cudaMemcpyDeviceToHost);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
            err = cudaFree(device_allocations[i]);
        if (err != cudaSuccess) {
            fprintf(stderr, "CUDA Error @ %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err));
            exit(3);
        }
        }
    }

        __device__ void operator()(int chunk_idx) {
            for (int __dummy_iter = 0; __dummy_iter < 1; __dummy_iter++) {
                {
            int i, j, k, i_global, j_global;
            float temp_top[BS*BS] __attribute__ ((aligned (64)));
            float temp_left[BS*BS] __attribute__ ((aligned (64)));
            float sum[BS] __attribute__ ((aligned (64))) = {0.f};
            
            i_global = offset + BS * (1 +  chunk_idx/chunks_in_inter_row);
            j_global = offset + BS * (1 + chunk_idx%chunks_in_inter_row);

            for (i = 0; i < BS; i++) {
for (j =0; j < BS; j++){
                    temp_top[i*BS + j]  = a[size*(i + offset) + j + j_global ];
                    temp_left[i*BS + j] = a[size*(i + i_global) + offset + j];
                }
            }

            for (i = 0; i < BS; i++)
            {
                for (k=0; k < BS; k++) {
for (j = 0; j < BS; j++) {
                        sum[j] += temp_left[BS*i + k] * temp_top[BS*k + j];
                    }
                }
for (j = 0; j < BS; j++) {
                    BB((i+i_global),(j+j_global)) -= sum[j];
                    sum[j] = 0.f;
                }
            }
        }
            }
        }
};

void lud_omp(float *a, int size)
{
{
    int offset, chunk_idx, size_inter, chunks_in_inter_row, chunks_per_inter;

    for (offset = 0; offset < size - BS ; offset += BS)
    {
        // lu factorization of left-top corner block diagonal matrix 
        //
        lud_diagonal_omp(a, size, offset);
            
        size_inter = size - offset -  BS;
        chunks_in_inter_row  = size_inter/BS;
        
        // calculate perimeter block matrices
        // 
 { const int niters = (chunks_in_inter_row) - (0);
const int iters_offset = (0);
kernel_launcher("pragma52_omp_parallel", iters_offset, niters, pragma52_omp_parallel_hclib_async(a, size, offset, chunk_idx));
 } 
        
        // update interior block matrices
        //
        chunks_per_inter = chunks_in_inter_row*chunks_in_inter_row;

 { const int niters = (chunks_per_inter) - (0);
const int iters_offset = (0);
kernel_launcher("pragma106_omp_parallel", iters_offset, niters, pragma106_omp_parallel_hclib_async(offset, chunk_idx, chunks_in_inter_row, a, size));
 } 
    }

    lud_diagonal_omp(a, size, offset);
    }
} 
