#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "blas.h"

#ifndef MIN
#define MIN(x,y)  (  ((x) < (y))? (x) : (y) )
#endif

#ifndef MAX
#define MAX(x,y)  (  ((x) > (y))? (x) : (y) )
#endif

#ifndef ABS
#define ABS(x) (( (x) > 0 )? (x) : (-(x)) )
#endif

#ifndef MOD
#define MOD(x,y)  ((x) % (y))
#endif

#ifndef FALSE
#define FALSE (1 == 0)
#endif

#ifndef TRUE
#define TRUE (1 == 1)
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern 
void estimate_kron_cost( const int nrow_A,
                         const int ncol_A,
                         const int nnz_A, 
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B,
                         double *p_kron_nnz, 
                         double *p_kron_flops, 
                         int *p_imethod );

extern
void csr2den(  const int nrow_A,
               const int ncol_A, 
               const int arowptr[],
               const int acol[],
               const double aval[],

               double a_[] );

extern
void csr_den_kron_mult_method( 
                    const int imethod,
                    const char transA,
                    const char transB,

                    const int nrow_A,
                    const int ncol_A, 
                    const int arowptr[], 
                    const int acol[], 
                    const double aval[],

                    const int nrow_B,
                    const int ncol_B, 
                    const double b_[],

                    const double yin[], 
                          double xout[] );

extern
int csr_nnz( const int nrow_A,  
             const int arowptr[] );

extern
void csr_transpose( 
                   const int nrow_A, 
                   const int ncol_A, 
                   const int arowptr[], 
                   const int acol[], 
                   const double aval[],  

                         int atrowptr[], 
                         int atcol[], 
                         double atval[] );

extern
void csr_kron_mult_method( 
                    const int imethod,
                    const char transA,
                    const char transB,
                   
                    const int nrow_A,
                    const int ncol_A, 
                    const int arowptr[], 
                    const int acol[], 
                    const double aval[],

                    const int nrow_B,
                    const int ncol_B, 
                    const int browptr[], 
                    const int bcol[], 
                    const double bval[],

                    const double yin[], 
                          double xout[] );



extern
void csr_matmul_post(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const int arowptr[],
                     const int acol[],
                     const double aval[],

                     const int nrow_Y,
                     const int ncol_Y,
                     const double yin[],

                     const int nrow_X,
                     const int ncol_X,
                     double xout[] );

extern
void csr_matmul_pre( const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const int arowptr[],
                     const int acol[], 
                     const double aval[],

                     const int nrow_Y,
                     const int ncol_Y,
                     const double yin[],

                     const int nrow_X,
                     const int ncol_X,
                     double xout[] );

extern
void csr_submatrix( const int nrow_A,
                    const int ncol_A,
                    const int arowptr[],
                    const int acol[],
                    const double aval[],

                    const int nrow_B,
                    const int ncol_B,
                    const int max_nnz,

                    const int rindex[],
                    const int cindex[],

                    int browptr[],
                    int bcol[],
                    double bval[] );


extern
void csr_kron_submatrix( 
         const int nrow_A,
         const int ncol_A,
         const int arowptr[],
         const int acol[],
         const double aval[],

         const int nrow_B,
         const int ncol_B,
         const int browptr[],
         const int bcol[],
         const double bval[],
         
         const int nrindex, 
         const int ncindex, 
         const int max_nnz,
         const int rindex[],
         const int cindex[],

         int hrowptr[],
         int hcol[],
         double hval[] );


extern
int csc_nnz( const int ncol_A,
             const int acolptr[] );

extern
void csc_matmul_pre( const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const int acolptr[],
                     const int arow[],
                     const double aval[],

                     const int nrow_Y,
                     const int ncol_Y,
                     const double yin[],

                     const int nrow_X,
                     const int ncol_X,
                     double xout[] );


extern
void csc_matmul_post(const char trans_A, 
                     const int nrow_A,
                     const int ncol_A,
                     const int acolptr[],
                     const int arow[],
                     const double aval[],

                     const int nrow_Y, 
                     const int ncol_Y,
                     const double yin[],

                     const int nrow_X,
                     const int ncol_X,
                     double xout[] );


extern
void csc_kron_mult_method( 
                    const int imethod,
                    const int nrow_A,
                    const int ncol_A, 
                    const int acolptr[], 
                    const int arow[], 
                    const double aval[],

                    const int nrow_B,
                    const int ncol_B, 
                    const int bcolptr[], 
                    const int brow[], 
                    const double bval[],

                    const double yin[], 
                          double xout[] );


extern
void csc_kron_mult( const int nrow_A,
                    const int ncol_A,
                    const int acolptr[],
                    const int arow[],
                    const double aval[],

                    const int nrow_B,
                    const int ncol_B,
                    const int bcolptr[],
                    const int brow[],
                    const double bval[],

                    const double yin[],
                          double xout[] );


extern
void coord2csr( 
              const int nrow_A, 
              const int ncol_A, 
              const int nnz, 
              const int ilist[], 
              const int jlist[], 
              const double alist[],
              int arowptr[], 
              int acol[], 
              double aval[] );



extern
void den_csr_kron_mult_method( 
                    const int imethod,
                    const char transA,
                    const char transB,

                    const int nrow_A,
                    const int ncol_A, 
                    const double a_[],

                    const int nrow_B,
                    const int ncol_B, 
                    const int browptr[], 
                    const int bcol[], 
                    const double bval[],

                    const double yin[], 
                          double xout[] );

extern
void den_copymat( const int nrow, 
                  const int ncol, 
                  const int asrc_[], 
                        int bdest_[] );

extern
void den_zeros( const int nrow_A,
                const int ncol_A,
                      double a_[] );

extern
void den_transpose( const int nrow_A, 
                    const int ncol_A, 
                    const double a_[],  
                          double at_[] );

extern
void den_gen_matrix( const int nrow_A, 
                     const int ncol_A, 
                     const double threshold, 
                           double a_[] );

extern
void den_matmul_pre( const char trans_A, 
                     const int nrow_A,
                     const int ncol_A, 
                     const double a_[],

                     const int nrow_Y, 
                     const int ncol_Y, 
                     const double yin[],

                     const int nrow_X, 
                     const int ncol_X, 
                     double xout[] );



extern
void den_matmul_post( 
                     const char trans_A, 
                     const int nrow_A,
                     const int ncol_A, 
                     const double a_[],

                     const int nrow_Y, 
                     const int ncol_Y, 
                     const double yin[],

                     const int nrow_X, 
                     const int ncol_X, 
                     double xout[] );



extern
void den_kron_submatrix(
        const int nrow_A,
        const int ncol_A,
        const double a_[],

        const int nrow_B,
        const int ncol_B,
        const double b_[],

        const int nrindex, 
        const int ncindex,
        const int rindex[], 
        const int cindex[],

        double c_[] );

extern
void den_kron_mult_method( 
                    const int imethod,
                    const char transA,
                    const char transB,

                    const int nrow_A,
                    const int ncol_A, 
                    const double a_[],

                    const int nrow_B,
                    const int ncol_B, 
                    const double b_[],

                    const double yin[], 
                          double xout[] );

extern
int den_nnz( const int nrow_A,
             const int ncol_A,
             const double  a_[] );


extern
void den_kron_form( const int nrow_A,
                    const int ncol_A, 
                    const double a_[],

                    const int nrow_B,
                    const int ncol_B, 
                    const double b_[],

                          double c_[] );


extern
void den_submatrix( const int nrow_A, 
                    const int ncol_A, 
                    const double a_[],

                    const int nrindex, 
                    const int ncindex,
                    const int rindex[],  
                    const int cindex[],

                    double c_[] );

extern
void den2csr( const int nrow_A, 
              const int ncol_A, 
              const double a_[], 

              const int max_nnz, 
                    int arowptr[], 
                    int acol[], 
                    double aval[] );





#ifdef __cplusplus
}
#endif

#endif

