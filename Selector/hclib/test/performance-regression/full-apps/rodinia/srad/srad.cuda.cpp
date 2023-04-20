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
// srad.cpp : Defines the entry point for the console application.
//

//#define OUTPUT


#define OPEN
#define	ITERATION
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

void random_matrix(float *I, int rows, int cols);

void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <rows> <cols> <y1> <y2> <x1> <x2> <no. of threads><lamda> <no. of iter>\n", argv[0]);
	fprintf(stderr, "\t<rows>   - number of rows\n");
	fprintf(stderr, "\t<cols>    - number of cols\n");
	fprintf(stderr, "\t<y1> 	 - y1 value of the speckle\n");
	fprintf(stderr, "\t<y2>      - y2 value of the speckle\n");
	fprintf(stderr, "\t<x1>       - x1 value of the speckle\n");
	fprintf(stderr, "\t<x2>       - x2 value of the speckle\n");
	fprintf(stderr, "\t<no. of threads>  - no. of threads\n");
	fprintf(stderr, "\t<lamda>   - lambda (0,1)\n");
	fprintf(stderr, "\t<no. of iter>   - number of iterations\n");
	
	exit(1);
}

class pragma125_omp_parallel_hclib_async {
    private:
        void **host_allocations;
        size_t *host_allocation_sizes;
        unsigned nallocations;
        void **device_allocations;
    volatile int cols;
    int k;
    float Jc;
    float* volatile J;
    float* volatile h_J;
    float* volatile dN;
    float* volatile h_dN;
    int* volatile iN;
    int* volatile h_iN;
    float* volatile dS;
    float* volatile h_dS;
    int* volatile iS;
    int* volatile h_iS;
    float* volatile dW;
    float* volatile h_dW;
    int* volatile jW;
    int* volatile h_jW;
    float* volatile dE;
    float* volatile h_dE;
    int* volatile jE;
    int* volatile h_jE;
    float G2;
    float L;
    float num;
    float den;
    float qsqr;
    volatile float q0sqr;
    float* volatile c;
    float* volatile h_c;

    public:
        pragma125_omp_parallel_hclib_async(int set_cols,
                int set_k,
                float set_Jc,
                float* set_J,
                float* set_dN,
                int* set_iN,
                float* set_dS,
                int* set_iS,
                float* set_dW,
                int* set_jW,
                float* set_dE,
                int* set_jE,
                float set_G2,
                float set_L,
                float set_num,
                float set_den,
                float set_qsqr,
                float set_q0sqr,
                float* set_c) {
            cols = set_cols;
            k = set_k;
            Jc = set_Jc;
            h_J = set_J;
            h_dN = set_dN;
            h_iN = set_iN;
            h_dS = set_dS;
            h_iS = set_iS;
            h_dW = set_dW;
            h_jW = set_jW;
            h_dE = set_dE;
            h_jE = set_jE;
            G2 = set_G2;
            L = set_L;
            num = set_num;
            den = set_den;
            qsqr = set_qsqr;
            q0sqr = set_q0sqr;
            h_c = set_c;

        }

    void transfer_to_device() {
        int i;
        cudaError_t err;

        J = NULL;
        dN = NULL;
        iN = NULL;
        dS = NULL;
        iS = NULL;
        dW = NULL;
        jW = NULL;
        dE = NULL;
        jE = NULL;
        c = NULL;

        get_underlying_allocations(&host_allocations, &host_allocation_sizes, &nallocations, 10, h_J, h_dN, h_iN, h_dS, h_iS, h_dW, h_jW, h_dE, h_jE, h_c);
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
            if (J == NULL && (char *)h_J >= (char *)host_allocations[i] && ((char *)h_J - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_J - (char *)host_allocations[i]);
                memcpy((void *)(&J), (void *)(&tmp), sizeof(void *));
            }
            if (dN == NULL && (char *)h_dN >= (char *)host_allocations[i] && ((char *)h_dN - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dN - (char *)host_allocations[i]);
                memcpy((void *)(&dN), (void *)(&tmp), sizeof(void *));
            }
            if (iN == NULL && (char *)h_iN >= (char *)host_allocations[i] && ((char *)h_iN - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_iN - (char *)host_allocations[i]);
                memcpy((void *)(&iN), (void *)(&tmp), sizeof(void *));
            }
            if (dS == NULL && (char *)h_dS >= (char *)host_allocations[i] && ((char *)h_dS - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dS - (char *)host_allocations[i]);
                memcpy((void *)(&dS), (void *)(&tmp), sizeof(void *));
            }
            if (iS == NULL && (char *)h_iS >= (char *)host_allocations[i] && ((char *)h_iS - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_iS - (char *)host_allocations[i]);
                memcpy((void *)(&iS), (void *)(&tmp), sizeof(void *));
            }
            if (dW == NULL && (char *)h_dW >= (char *)host_allocations[i] && ((char *)h_dW - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dW - (char *)host_allocations[i]);
                memcpy((void *)(&dW), (void *)(&tmp), sizeof(void *));
            }
            if (jW == NULL && (char *)h_jW >= (char *)host_allocations[i] && ((char *)h_jW - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_jW - (char *)host_allocations[i]);
                memcpy((void *)(&jW), (void *)(&tmp), sizeof(void *));
            }
            if (dE == NULL && (char *)h_dE >= (char *)host_allocations[i] && ((char *)h_dE - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dE - (char *)host_allocations[i]);
                memcpy((void *)(&dE), (void *)(&tmp), sizeof(void *));
            }
            if (jE == NULL && (char *)h_jE >= (char *)host_allocations[i] && ((char *)h_jE - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_jE - (char *)host_allocations[i]);
                memcpy((void *)(&jE), (void *)(&tmp), sizeof(void *));
            }
            if (c == NULL && (char *)h_c >= (char *)host_allocations[i] && ((char *)h_c - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_c - (char *)host_allocations[i]);
                memcpy((void *)(&c), (void *)(&tmp), sizeof(void *));
            }
        }

        assert(J || h_J == NULL);
        assert(dN || h_dN == NULL);
        assert(iN || h_iN == NULL);
        assert(dS || h_dS == NULL);
        assert(iS || h_iS == NULL);
        assert(dW || h_dW == NULL);
        assert(jW || h_jW == NULL);
        assert(dE || h_dE == NULL);
        assert(jE || h_jE == NULL);
        assert(c || h_c == NULL);

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

        __device__ void operator()(int i) {
            for (int __dummy_iter = 0; __dummy_iter < 1; __dummy_iter++) {
                {
            for (int j = 0; j < cols; j++) { 
		
				k = i * cols + j;
				Jc = J[k];
 
				// directional derivates
                dN[k] = J[iN[i] * cols + j] - Jc;
                dS[k] = J[iS[i] * cols + j] - Jc;
                dW[k] = J[i * cols + jW[j]] - Jc;
                dE[k] = J[i * cols + jE[j]] - Jc;
			
                G2 = (dN[k]*dN[k] + dS[k]*dS[k] 
                    + dW[k]*dW[k] + dE[k]*dE[k]) / (Jc*Jc);

   		        L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;

				num  = (0.5*G2) - ((1.0/16.0)*(L*L)) ;
                den  = 1 + (.25*L);
                qsqr = num/(den*den);
 
                // diffusion coefficent (equ 33)
                den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;
                c[k] = 1.0 / (1.0+den) ;
                
                // saturate diffusion coefficent
                if (c[k] < 0) {c[k] = 0;}
                else if (c[k] > 1) {c[k] = 1;}
   
		}
  
    }
            }
        }
};

class pragma158_omp_parallel_hclib_async {
    private:
        void **host_allocations;
        size_t *host_allocation_sizes;
        unsigned nallocations;
        void **device_allocations;
    volatile int cols;
    int k;
    float cN;
    float* volatile c;
    float* volatile h_c;
    float cS;
    int* volatile iS;
    int* volatile h_iS;
    float cW;
    float cE;
    int* volatile jE;
    int* volatile h_jE;
    float D;
    float* volatile dN;
    float* volatile h_dN;
    float* volatile dS;
    float* volatile h_dS;
    float* volatile dW;
    float* volatile h_dW;
    float* volatile dE;
    float* volatile h_dE;
    float* volatile J;
    float* volatile h_J;
    volatile float lambda;

    public:
        pragma158_omp_parallel_hclib_async(int set_cols,
                int set_k,
                float set_cN,
                float* set_c,
                float set_cS,
                int* set_iS,
                float set_cW,
                float set_cE,
                int* set_jE,
                float set_D,
                float* set_dN,
                float* set_dS,
                float* set_dW,
                float* set_dE,
                float* set_J,
                float set_lambda) {
            cols = set_cols;
            k = set_k;
            cN = set_cN;
            h_c = set_c;
            cS = set_cS;
            h_iS = set_iS;
            cW = set_cW;
            cE = set_cE;
            h_jE = set_jE;
            D = set_D;
            h_dN = set_dN;
            h_dS = set_dS;
            h_dW = set_dW;
            h_dE = set_dE;
            h_J = set_J;
            lambda = set_lambda;

        }

    void transfer_to_device() {
        int i;
        cudaError_t err;

        c = NULL;
        iS = NULL;
        jE = NULL;
        dN = NULL;
        dS = NULL;
        dW = NULL;
        dE = NULL;
        J = NULL;

        get_underlying_allocations(&host_allocations, &host_allocation_sizes, &nallocations, 8, h_c, h_iS, h_jE, h_dN, h_dS, h_dW, h_dE, h_J);
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
            if (c == NULL && (char *)h_c >= (char *)host_allocations[i] && ((char *)h_c - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_c - (char *)host_allocations[i]);
                memcpy((void *)(&c), (void *)(&tmp), sizeof(void *));
            }
            if (iS == NULL && (char *)h_iS >= (char *)host_allocations[i] && ((char *)h_iS - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_iS - (char *)host_allocations[i]);
                memcpy((void *)(&iS), (void *)(&tmp), sizeof(void *));
            }
            if (jE == NULL && (char *)h_jE >= (char *)host_allocations[i] && ((char *)h_jE - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_jE - (char *)host_allocations[i]);
                memcpy((void *)(&jE), (void *)(&tmp), sizeof(void *));
            }
            if (dN == NULL && (char *)h_dN >= (char *)host_allocations[i] && ((char *)h_dN - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dN - (char *)host_allocations[i]);
                memcpy((void *)(&dN), (void *)(&tmp), sizeof(void *));
            }
            if (dS == NULL && (char *)h_dS >= (char *)host_allocations[i] && ((char *)h_dS - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dS - (char *)host_allocations[i]);
                memcpy((void *)(&dS), (void *)(&tmp), sizeof(void *));
            }
            if (dW == NULL && (char *)h_dW >= (char *)host_allocations[i] && ((char *)h_dW - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dW - (char *)host_allocations[i]);
                memcpy((void *)(&dW), (void *)(&tmp), sizeof(void *));
            }
            if (dE == NULL && (char *)h_dE >= (char *)host_allocations[i] && ((char *)h_dE - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_dE - (char *)host_allocations[i]);
                memcpy((void *)(&dE), (void *)(&tmp), sizeof(void *));
            }
            if (J == NULL && (char *)h_J >= (char *)host_allocations[i] && ((char *)h_J - (char *)host_allocations[i]) < host_allocation_sizes[i]) {
                char *tmp = (char *)device_allocations[i] + ((char *)h_J - (char *)host_allocations[i]);
                memcpy((void *)(&J), (void *)(&tmp), sizeof(void *));
            }
        }

        assert(c || h_c == NULL);
        assert(iS || h_iS == NULL);
        assert(jE || h_jE == NULL);
        assert(dN || h_dN == NULL);
        assert(dS || h_dS == NULL);
        assert(dW || h_dW == NULL);
        assert(dE || h_dE == NULL);
        assert(J || h_J == NULL);

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

        __device__ void operator()(int i) {
            for (int __dummy_iter = 0; __dummy_iter < 1; __dummy_iter++) {
                {
            for (int j = 0; j < cols; j++) {        

                // current index
                k = i * cols + j;
                
                // diffusion coefficent
					cN = c[k];
					cS = c[iS[i] * cols + j];
					cW = c[k];
					cE = c[i * cols + jE[j]];

                // divergence (equ 58)
                D = cN * dN[k] + cS * dS[k] + cW * dW[k] + cE * dE[k];
                
                // image update (equ 61)
                J[k] = J[k] + 0.25*lambda*D;
                #ifdef OUTPUT
                //printf("%.5f ", J[k]); 
                #endif //output
            }
	            #ifdef OUTPUT
                //printf("\n"); 
                #endif //output
	     }
            }
        }
};

int main(int argc, char* argv[])
{   
	int rows, cols, size_I, size_R, niter = 10, iter, k;
    float *I, *J, q0sqr, sum, sum2, tmp, meanROI,varROI ;
	float Jc, G2, L, num, den, qsqr;
	int *iN,*iS,*jE,*jW;
	float *dN,*dS,*dW,*dE;
	int r1, r2, c1, c2;
	float cN,cS,cW,cE;
	float *c, D;
	float lambda;
	int i, j;
    int nthreads;

	if (argc == 10)
	{
		rows = atoi(argv[1]); //number of rows in the domain
		cols = atoi(argv[2]); //number of cols in the domain
		if ((rows%16!=0) || (cols%16!=0)){
			fprintf(stderr, "rows and cols must be multiples of 16\n");
			exit(1);
		}
		r1   = atoi(argv[3]); //y1 position of the speckle
		r2   = atoi(argv[4]); //y2 position of the speckle
		c1   = atoi(argv[5]); //x1 position of the speckle
		c2   = atoi(argv[6]); //x2 position of the speckle
		nthreads = atoi(argv[7]); // number of threads
		lambda = atof(argv[8]); //Lambda value
		niter = atoi(argv[9]); //number of iterations
	}
    else{
		usage(argc, argv);
    }


	size_I = cols * rows;
    size_R = (r2-r1+1)*(c2-c1+1);   

	I = (float *)malloc( size_I * sizeof(float) );
    J = (float *)malloc( size_I * sizeof(float) );
	c  = (float *)malloc(sizeof(float)* size_I) ;

    iN = (int *)malloc(sizeof(unsigned int*) * rows) ;
    iS = (int *)malloc(sizeof(unsigned int*) * rows) ;
    jW = (int *)malloc(sizeof(unsigned int*) * cols) ;
    jE = (int *)malloc(sizeof(unsigned int*) * cols) ;    


	dN = (float *)malloc(sizeof(float)* size_I) ;
    dS = (float *)malloc(sizeof(float)* size_I) ;
    dW = (float *)malloc(sizeof(float)* size_I) ;
    dE = (float *)malloc(sizeof(float)* size_I) ;    
    

    for (int i=0; i< rows; i++) {
        iN[i] = i-1;
        iS[i] = i+1;
    }    
    for (int j=0; j< cols; j++) {
        jW[j] = j-1;
        jE[j] = j+1;
    }
    iN[0]    = 0;
    iS[rows-1] = rows-1;
    jW[0]    = 0;
    jE[cols-1] = cols-1;
	
	printf("Randomizing the input matrix\n");

    random_matrix(I, rows, cols);

    for (k = 0;  k < size_I; k++ ) {
     	J[k] = (float)exp(I[k]) ;
    }
   
	printf("Start the SRAD main loop\n");

for (iter=0; iter< niter; iter++){
		sum=0; sum2=0;     
		for (i=r1; i<=r2; i++) {
            for (j=c1; j<=c2; j++) {
                tmp   = J[i * cols + j];
                sum  += tmp ;
                sum2 += tmp*tmp;
            }
        }
        meanROI = sum / size_R;
        varROI  = (sum2 / size_R) - meanROI*meanROI;
        q0sqr   = varROI / (meanROI*meanROI);
		

 { const int niters = (rows) - (0);
const int iters_offset = (0);
kernel_launcher("pragma125_omp_parallel", iters_offset, niters, pragma125_omp_parallel_hclib_async(cols, k, Jc, J, dN, iN, dS, iS, dW, jW, dE, jE, G2, L, num, den, qsqr, q0sqr, c));
 } 
 { const int niters = (rows) - (0);
const int iters_offset = (0);
kernel_launcher("pragma158_omp_parallel", iters_offset, niters, pragma158_omp_parallel_hclib_async(cols, k, cN, c, cS, iS, cW, cE, jE, D, dN, dS, dW, dE, J, lambda));
 } 

	}


#ifdef OUTPUT
	  for( int i = 0 ; i < rows ; i++){
		for ( int j = 0 ; j < cols ; j++){

         printf("%.5f ", J[i * cols + j]); 
    
		}
         printf("\n"); 
   }
#endif 

	printf("Computation Done\n");

	free(I);
	free(J);
	free(iN); free(iS); free(jW); free(jE);
    free(dN); free(dS); free(dW); free(dE);

	free(c);
	return 0;
} 




void random_matrix(float *I, int rows, int cols){

	srand(7);
	
	for( int i = 0 ; i < rows ; i++){
		for ( int j = 0 ; j < cols ; j++){
		 I[i * cols + j] = rand()/(float)RAND_MAX ;
		 #ifdef OUTPUT
         //printf("%g ", I[i * cols + j]); 
         #endif 
		}
		 #ifdef OUTPUT
         //printf("\n"); 
         #endif 
	}

}

