#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "Matrix.h"
#include "MatrixNonOwned.h"
#include "CrsMatrix.h"

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
void csr_den_kron_mult_method(const int imethod,
                    const char transA,
                    const char transB,

                    const PsimagLite::CrsMatrix<double>& a_,
                    const PsimagLite::Matrix<double>& b_,

                    const double* yin,
                          double* xout );

extern
int csr_nnz(const PsimagLite::CrsMatrix<double>&);

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
void csr_kron_mult_method(const int imethod,
                    const char transA,
                    const char transB,

                    const PsimagLite::CrsMatrix<double>& a,
                    const PsimagLite::CrsMatrix<double>& b,

                    const PsimagLite::MatrixNonOwned<const double>& yin,
                          PsimagLite::MatrixNonOwned<double>& xout);



extern
void csr_matmul_post(const char trans_A,

                     const PsimagLite::CrsMatrix<double>&,

                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);

extern
void csr_matmul_pre( const char trans_A,

                     const PsimagLite::CrsMatrix<double>&,

                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);

extern
void csr_submatrix(const PsimagLite::CrsMatrix<double>& a,
                   const int nrow_B,
				           const int ncol_B,
                    const int max_nnz,

                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,
                    PsimagLite::CrsMatrix<double>& b);


extern
void csr_kron_submatrix(const PsimagLite::CrsMatrix<double>& a,
         const PsimagLite::CrsMatrix<double>& b,

         const int nrindex,
         const int ncindex,
         const int max_nnz,
         const PsimagLite::Vector<int>::Type& rindex,
         const PsimagLite::Vector<int>::Type& cindex,
         PsimagLite::CrsMatrix<double>& h);


extern
int csc_nnz(const int ncol_A,
             const PsimagLite::Vector<int>::Type& acolptr );

extern
void csc_matmul_pre( const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<double>::Type& aval,

                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<double>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<double>& xout );


extern
void csc_matmul_post(const char trans_A, 
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<double>::Type& aval,

                     const int nrow_Y, 
                     const int ncol_Y,
                     const PsimagLite::Matrix<double>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<double>& xout );


extern
void csc_kron_mult_method(const int imethod,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Vector<int>::Type& acolptr,
                    const PsimagLite::Vector<int>::Type& arow,
                    const PsimagLite::Vector<double>::Type& aval,

                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Vector<int>::Type& bcolptr,
                    const PsimagLite::Vector<int>::Type& brow,
                    const PsimagLite::Vector<double>::Type& bval,

                    const PsimagLite::Matrix<double>& yin,
                          PsimagLite::Matrix<double>& xout );


extern
void csc_kron_mult(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Vector<int>::Type& acolptr,
                    const PsimagLite::Vector<int>::Type& arow,
                    const PsimagLite::Vector<double>::Type& aval,

                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Vector<int>::Type& bcolptr,
                    const PsimagLite::Vector<int>::Type& brow,
                    const PsimagLite::Vector<double>::Type& bval,

                    const PsimagLite::Matrix<double>& yin,
                          PsimagLite::Matrix<double>& xout );


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
void den_csr_kron_mult_method(const int imethod,
                    const char transA,
                    const char transB,
                    const PsimagLite::Matrix<double>& a_,
                    const PsimagLite::CrsMatrix<double>& b,

                    const double* yin,
                          double* xout_);

extern
void den_copymat( const int nrow, 
                  const int ncol, 
                  const int asrc_[], 
                        int bdest_[] );

extern
void den_zeros( const int nrow_A,
                const int ncol_A,
                      PsimagLite::Matrix<double>& a_);

extern
void den_transpose( const int nrow_A, 
                    const int ncol_A, 
                    const double a_[],  
                          double at_[] );

extern
void den_gen_matrix( const int nrow_A, 
                     const int ncol_A, 
                     const double threshold, 
                           PsimagLite::Matrix<double>& a_);

extern
void den_matmul_pre(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<double>& a_,

                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);



extern
void den_matmul_post(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<double>& a_,

                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);



extern
void den_kron_submatrix(
        const int nrow_A,
        const int ncol_A,
        const PsimagLite::Matrix<double>& a_,

        const int nrow_B,
        const int ncol_B,
        const PsimagLite::Matrix<double>& b_,

        const int nrindex, 
        const int ncindex,
        const PsimagLite::Vector<int>::Type& rindex,
        const PsimagLite::Vector<int>::Type& cindex,

        PsimagLite::Matrix<double>& c_ );

extern
void den_kron_mult_method(const int imethod,
                    const char transA,
                    const char transB,

                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<double>& a_,

                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<double>& b_,

                    const double* yin_,
                          double* xout_ );

extern
int den_nnz(const PsimagLite::Matrix<double>&);


extern
void den_kron_form( const int nrow_A,
                    const int ncol_A, 
                    const PsimagLite::Matrix<double>& a_,

                    const int nrow_B,
                    const int ncol_B, 
                    const PsimagLite::Matrix<double>& b_,

                          PsimagLite::Matrix<double>& c_);


extern
void den_submatrix( const int nrow_A, 
                    const int ncol_A, 
                    const PsimagLite::Matrix<double>& a_,

                    const int nrindex, 
                    const int ncindex,
                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,

                    PsimagLite::Matrix<double>& c_);


#ifdef __cplusplus
}
#endif

#endif

