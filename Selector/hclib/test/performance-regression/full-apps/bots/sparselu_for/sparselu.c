#include "hclib.h"
#ifdef __cplusplus
#include "hclib_cpp.h"
#include "hclib_system.h"
#ifdef __CUDACC__
#include "hclib_cuda.h"
#endif
#endif
/**********************************************************************************************/
/*  This program is part of the Barcelona OpenMP Tasks Suite                                  */
/*  Copyright (C) 2009 Barcelona Supercomputing Center - Centro Nacional de Supercomputacion  */
/*  Copyright (C) 2009 Universitat Politecnica de Catalunya                                   */
/*                                                                                            */
/*  This program is free software; you can redistribute it and/or modify                      */
/*  it under the terms of the GNU General Public License as published by                      */
/*  the Free Software Foundation; either version 2 of the License, or                         */
/*  (at your option) any later version.                                                       */
/*                                                                                            */
/*  This program is distributed in the hope that it will be useful,                           */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of                            */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             */
/*  GNU General Public License for more details.                                              */
/*                                                                                            */
/*  You should have received a copy of the GNU General Public License                         */
/*  along with this program; if not, write to the Free Software                               */
/*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA            */
/**********************************************************************************************/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include <libgen.h>
#include "bots.h"
#include "sparselu.h"

/***********************************************************************
 * checkmat: 
 **********************************************************************/
int checkmat (float *M, float *N)
{
   int i, j;
   float r_err;

   for (i = 0; i < bots_arg_size_1; i++) 
   {
      for (j = 0; j < bots_arg_size_1; j++) 
      {
         r_err = M[i*bots_arg_size_1+j] - N[i*bots_arg_size_1+j];
         if ( r_err == 0.0 ) continue;

         if (r_err < 0.0 ) r_err = -r_err;

         if ( M[i*bots_arg_size_1+j] == 0 )
         {
           bots_message("Checking failure: A[%d][%d]=%f  B[%d][%d]=%f; \n",
                    i,j, M[i*bots_arg_size_1+j], i,j, N[i*bots_arg_size_1+j]);
           return FALSE;
         }
         r_err = r_err / M[i*bots_arg_size_1+j];
         if(r_err > EPSILON)
         {
            bots_message("Checking failure: A[%d][%d]=%f  B[%d][%d]=%f; Relative Error=%f\n",
                    i,j, M[i*bots_arg_size_1+j], i,j, N[i*bots_arg_size_1+j], r_err);
            return FALSE;
         }
      }
   }
   return TRUE;
}
/***********************************************************************
 * genmat: 
 **********************************************************************/
void genmat (float *M[])
{
   int null_entry, init_val, i, j, ii, jj;
   float *p;
   int a=0,b=0;

   init_val = 1325;

   /* generating the structure */
   for (ii=0; ii < bots_arg_size; ii++)
   {
      for (jj=0; jj < bots_arg_size; jj++)
      {
         /* computing null entries */
         null_entry=FALSE;
         if ((ii<jj) && (ii%3 !=0)) null_entry = TRUE;
         if ((ii>jj) && (jj%3 !=0)) null_entry = TRUE;
	 if (ii%2==1) null_entry = TRUE;
	 if (jj%2==1) null_entry = TRUE;
	 if (ii==jj) null_entry = FALSE;
	 if (ii==jj-1) null_entry = FALSE;
         if (ii-1 == jj) null_entry = FALSE; 
         /* allocating matrix */
         if (null_entry == FALSE){
            a++;
            M[ii*bots_arg_size+jj] = (float *) malloc(bots_arg_size_1*bots_arg_size_1*sizeof(float));
	    if ((M[ii*bots_arg_size+jj] == NULL))
            {
               bots_message("Error: Out of memory\n");
               exit(101);
            }
            /* initializing matrix */
            p = M[ii*bots_arg_size+jj];
            for (i = 0; i < bots_arg_size_1; i++) 
            {
               for (j = 0; j < bots_arg_size_1; j++)
               {
	            init_val = (3125 * init_val) % 65536;
      	            (*p) = (float)((init_val - 32768.0) / 16384.0);
                    p++;
               }
            }
         }
         else
         {
            b++;
            M[ii*bots_arg_size+jj] = NULL;
         }
      }
   }
   bots_debug("allo = %d, no = %d, total = %d, factor = %f\n",a,b,a+b,(float)((float)a/(float)(a+b)));
}
/***********************************************************************
 * print_structure: 
 **********************************************************************/
void print_structure(char *name, float *M[])
{
   int ii, jj;
   bots_message("Structure for matrix %s @ 0x%p\n",name, M);
   for (ii = 0; ii < bots_arg_size; ii++) {
     for (jj = 0; jj < bots_arg_size; jj++) {
        if (M[ii*bots_arg_size+jj]!=NULL) {bots_message("x");}
        else bots_message(" ");
     }
     bots_message("\n");
   }
   bots_message("\n");
}
/***********************************************************************
 * allocate_clean_block: 
 **********************************************************************/
float * allocate_clean_block()
{
  int i,j;
  float *p, *q;

  p = (float *) malloc(bots_arg_size_1*bots_arg_size_1*sizeof(float));
  q=p;
  if (p!=NULL){
     for (i = 0; i < bots_arg_size_1; i++) 
        for (j = 0; j < bots_arg_size_1; j++){(*p)=0.0; p++;}
	
  }
  else
  {
      bots_message("Error: Out of memory\n");
      exit (101);
  }
  return (q);
}

/***********************************************************************
 * lu0: 
 **********************************************************************/
void lu0(float *diag)
{
   int i, j, k;

   for (k=0; k<bots_arg_size_1; k++)
      for (i=k+1; i<bots_arg_size_1; i++)
      {
         diag[i*bots_arg_size_1+k] = diag[i*bots_arg_size_1+k] / diag[k*bots_arg_size_1+k];
         for (j=k+1; j<bots_arg_size_1; j++)
            diag[i*bots_arg_size_1+j] = diag[i*bots_arg_size_1+j] - diag[i*bots_arg_size_1+k] * diag[k*bots_arg_size_1+j];
      }
}

/***********************************************************************
 * bdiv: 
 **********************************************************************/
void bdiv(float *diag, float *row)
{
   int i, j, k;
   for (i=0; i<bots_arg_size_1; i++)
      for (k=0; k<bots_arg_size_1; k++)
      {
         row[i*bots_arg_size_1+k] = row[i*bots_arg_size_1+k] / diag[k*bots_arg_size_1+k];
         for (j=k+1; j<bots_arg_size_1; j++)
            row[i*bots_arg_size_1+j] = row[i*bots_arg_size_1+j] - row[i*bots_arg_size_1+k]*diag[k*bots_arg_size_1+j];
      }
}
/***********************************************************************
 * bmod: 
 **********************************************************************/
void bmod(float *row, float *col, float *inner)
{
   int i, j, k;
   for (i=0; i<bots_arg_size_1; i++)
      for (j=0; j<bots_arg_size_1; j++)
         for (k=0; k<bots_arg_size_1; k++)
            inner[i*bots_arg_size_1+j] = inner[i*bots_arg_size_1+j] - row[i*bots_arg_size_1+k]*col[k*bots_arg_size_1+j];
}
/***********************************************************************
 * fwd: 
 **********************************************************************/
void fwd(float *diag, float *col)
{
   int i, j, k;
   for (j=0; j<bots_arg_size_1; j++)
      for (k=0; k<bots_arg_size_1; k++) 
         for (i=k+1; i<bots_arg_size_1; i++)
            col[i*bots_arg_size_1+j] = col[i*bots_arg_size_1+j] - diag[i*bots_arg_size_1+k]*col[k*bots_arg_size_1+j];
}

void sparselu_init (float ***pBENCH, char *pass)
{
   *pBENCH = (float **) malloc(bots_arg_size*bots_arg_size*sizeof(float *));
   genmat(*pBENCH);
   print_structure(pass, *pBENCH);
}


void sparselu_seq_call(float **BENCH)
{
   int ii, jj, kk;

   for (kk=0; kk<bots_arg_size; kk++)
   {
      lu0(BENCH[kk*bots_arg_size+kk]);
      for (jj=kk+1; jj<bots_arg_size; jj++)
         if (BENCH[kk*bots_arg_size+jj] != NULL)
         {
            fwd(BENCH[kk*bots_arg_size+kk], BENCH[kk*bots_arg_size+jj]);
         }
      for (ii=kk+1; ii<bots_arg_size; ii++) 
         if (BENCH[ii*bots_arg_size+kk] != NULL)
         {
            bdiv (BENCH[kk*bots_arg_size+kk], BENCH[ii*bots_arg_size+kk]);
         }
      for (ii=kk+1; ii<bots_arg_size; ii++)
         if (BENCH[ii*bots_arg_size+kk] != NULL)
            for (jj=kk+1; jj<bots_arg_size; jj++)
               if (BENCH[kk*bots_arg_size+jj] != NULL)
               {
                     if (BENCH[ii*bots_arg_size+jj]==NULL) BENCH[ii*bots_arg_size+jj] = allocate_clean_block();
                     bmod(BENCH[ii*bots_arg_size+kk], BENCH[kk*bots_arg_size+jj], BENCH[ii*bots_arg_size+jj]);
               }

   }
}

typedef struct _pragma266_omp_task {
    int (*ii_ptr);
    int jj;
    int kk;
    float (*(*(*BENCH_ptr)));
 } pragma266_omp_task;

typedef struct _pragma274_omp_task {
    int ii;
    int (*jj_ptr);
    int kk;
    float (*(*(*BENCH_ptr)));
 } pragma274_omp_task;

typedef struct _pragma285_omp_task {
    int ii;
    int jj;
    int kk;
    float (*(*(*BENCH_ptr)));
 } pragma285_omp_task;

static void pragma266_omp_task_hclib_async(void *____arg);
static void pragma274_omp_task_hclib_async(void *____arg);
static void pragma285_omp_task_hclib_async(void *____arg);
typedef struct _main_entrypoint_ctx {
    int ii;
    int jj;
    int kk;
    float (*(*BENCH));
 } main_entrypoint_ctx;


static void main_entrypoint(void *____arg) {
    main_entrypoint_ctx *ctx = (main_entrypoint_ctx *)____arg;
    int ii; ii = ctx->ii;
    int jj; jj = ctx->jj;
    int kk; kk = ctx->kk;
    float (*(*BENCH)); BENCH = ctx->BENCH;
{
   bots_message("Computing SparseLU Factorization (%dx%d matrix with %dx%d blocks) ",
           bots_arg_size,bots_arg_size,bots_arg_size_1,bots_arg_size_1);
hclib_start_finish(); for (kk=0; kk<bots_arg_size; kk++) 
   {
      lu0(BENCH[kk*bots_arg_size+kk]);

      for (jj=kk+1; jj<bots_arg_size; jj++)
         if (BENCH[kk*bots_arg_size+jj] != NULL)
         {
 { 
pragma266_omp_task *new_ctx = (pragma266_omp_task *)malloc(sizeof(pragma266_omp_task));
new_ctx->ii_ptr = &(ii);
new_ctx->jj = jj;
new_ctx->kk = kk;
new_ctx->BENCH_ptr = &(BENCH);
hclib_async(pragma266_omp_task_hclib_async, new_ctx, NO_FUTURE, ANY_PLACE);
 } 
         }
      for (ii=kk+1; ii<bots_arg_size; ii++) 
         if (BENCH[ii*bots_arg_size+kk] != NULL)
         {
 { 
pragma274_omp_task *new_ctx = (pragma274_omp_task *)malloc(sizeof(pragma274_omp_task));
new_ctx->ii = ii;
new_ctx->jj_ptr = &(jj);
new_ctx->kk = kk;
new_ctx->BENCH_ptr = &(BENCH);
hclib_async(pragma274_omp_task_hclib_async, new_ctx, NO_FUTURE, ANY_PLACE);
 } 
         }

      for (ii=kk+1; ii<bots_arg_size; ii++)
         if (BENCH[ii*bots_arg_size+kk] != NULL)
            for (jj=kk+1; jj<bots_arg_size; jj++)
               if (BENCH[kk*bots_arg_size+jj] != NULL)
               {
 { 
pragma285_omp_task *new_ctx = (pragma285_omp_task *)malloc(sizeof(pragma285_omp_task));
new_ctx->ii = ii;
new_ctx->jj = jj;
new_ctx->kk = kk;
new_ctx->BENCH_ptr = &(BENCH);
hclib_async(pragma285_omp_task_hclib_async, new_ctx, NO_FUTURE, ANY_PLACE);
 } 
               }

   } ; hclib_end_finish(); 
   bots_message(" completed!\n");
   } ;     free(____arg);
}

void sparselu_par_call(float **BENCH)
{
   int ii, jj, kk;
 
main_entrypoint_ctx *new_ctx = (main_entrypoint_ctx *)malloc(sizeof(main_entrypoint_ctx));
new_ctx->ii = ii;
new_ctx->jj = jj;
new_ctx->kk = kk;
new_ctx->BENCH = BENCH;
const char *deps[] = { "system" };
hclib_launch(main_entrypoint, new_ctx, deps, 1);

}  
static void pragma266_omp_task_hclib_async(void *____arg) {
    pragma266_omp_task *ctx = (pragma266_omp_task *)____arg;
    int jj; jj = ctx->jj;
    int kk; kk = ctx->kk;
    hclib_start_finish();
{
            fwd((*(ctx->BENCH_ptr))[kk*bots_arg_size+kk], (*(ctx->BENCH_ptr))[kk*bots_arg_size+jj]);
             } ;     ; hclib_end_finish_nonblocking();

    free(____arg);
}


static void pragma274_omp_task_hclib_async(void *____arg) {
    pragma274_omp_task *ctx = (pragma274_omp_task *)____arg;
    int ii; ii = ctx->ii;
    int kk; kk = ctx->kk;
    hclib_start_finish();
{
            bdiv ((*(ctx->BENCH_ptr))[kk*bots_arg_size+kk], (*(ctx->BENCH_ptr))[ii*bots_arg_size+kk]);
             } ;     ; hclib_end_finish_nonblocking();

    free(____arg);
}


static void pragma285_omp_task_hclib_async(void *____arg) {
    pragma285_omp_task *ctx = (pragma285_omp_task *)____arg;
    int ii; ii = ctx->ii;
    int jj; jj = ctx->jj;
    int kk; kk = ctx->kk;
    hclib_start_finish();
{
                     if ((*(ctx->BENCH_ptr))[ii*bots_arg_size+jj]==NULL) (*(ctx->BENCH_ptr))[ii*bots_arg_size+jj] = allocate_clean_block();
                     bmod((*(ctx->BENCH_ptr))[ii*bots_arg_size+kk], (*(ctx->BENCH_ptr))[kk*bots_arg_size+jj], (*(ctx->BENCH_ptr))[ii*bots_arg_size+jj]);
                   } ;     ; hclib_end_finish_nonblocking();

    free(____arg);
}



void sparselu_fini (float **BENCH, char *pass)
{
   print_structure(pass, BENCH);
}

int sparselu_check(float **SEQ, float **BENCH)
{
   int ii,jj,ok=1;

   for (ii=0; ((ii<bots_arg_size) && ok); ii++)
   {
      for (jj=0; ((jj<bots_arg_size) && ok); jj++)
      {
         if ((SEQ[ii*bots_arg_size+jj] == NULL) && (BENCH[ii*bots_arg_size+jj] != NULL)) ok = FALSE;
         if ((SEQ[ii*bots_arg_size+jj] != NULL) && (BENCH[ii*bots_arg_size+jj] == NULL)) ok = FALSE;
         if ((SEQ[ii*bots_arg_size+jj] != NULL) && (BENCH[ii*bots_arg_size+jj] != NULL))
            ok = checkmat(SEQ[ii*bots_arg_size+jj], BENCH[ii*bots_arg_size+jj]);
      }
   }
   if (ok) return BOTS_RESULT_SUCCESSFUL;
   else return BOTS_RESULT_UNSUCCESSFUL;
}

