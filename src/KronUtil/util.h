#ifndef UTIL_H
#define UTIL_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "KronUtil.h"
#include "MatrixNonOwned.h"

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

template<typename ComplexOrRealType>
void estimate_kron_cost( const int nrow_A,
                         const int ncol_A,
                         const int nnz_A,
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B,
                         ComplexOrRealType *p_kron_nnz,
                         ComplexOrRealType *p_kron_flops,
                         int *p_imethod );

template<typename ComplexOrRealType>
void csr_den_kron_mult_method(const int imethod,
                              const char transA,
                              const char transB,

                              const PsimagLite::CrsMatrix<ComplexOrRealType>& a_,
                              const PsimagLite::Matrix<ComplexOrRealType>& b_,

                              const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                              SizeType offsetY ,
                              typename PsimagLite::Vector<ComplexOrRealType>::Type& xout_,
                              SizeType offsetX);

template<typename ComplexOrRealType>
int csr_nnz(const PsimagLite::CrsMatrix<ComplexOrRealType>&);

template<typename ComplexOrRealType>
bool csr_is_eye(const PsimagLite::CrsMatrix<ComplexOrRealType>&);

template<typename ComplexOrRealType>
bool csr_is_zeros(const PsimagLite::CrsMatrix<ComplexOrRealType>&);

template<typename ComplexOrRealType>
void csr_transpose(const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const ComplexOrRealType aval[],
                   int atrowptr[],
                   int atcol[],
                   ComplexOrRealType atval[] );

template<typename ComplexOrRealType>
void csr_kron_mult_method(const int imethod,
                          const char transA,
                          const char transB,

                          const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                          const PsimagLite::CrsMatrix<ComplexOrRealType>& b,

                          const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                          PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout);



template<typename ComplexOrRealType>
void csr_matmul_post(const char trans_A,
                     const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout);

template<typename ComplexOrRealType>
void csr_matmul_pre( const char trans_A,

                     const PsimagLite::CrsMatrix<ComplexOrRealType>&,

                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,

                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout);

template<typename ComplexOrRealType>
void csr_submatrix(const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                   const int nrow_B,
                   const int ncol_B,
                   const int max_nnz,

                   const PsimagLite::Vector<int>::Type& rindex,
                   const PsimagLite::Vector<int>::Type& cindex,
                   PsimagLite::CrsMatrix<ComplexOrRealType>& b);


template<typename ComplexOrRealType>
void csr_eye(const int nrow_B,
             const int ncol_B,
             PsimagLite::CrsMatrix<ComplexOrRealType>& b);

template<typename ComplexOrRealType>
void csr_kron_submatrix(const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                        const PsimagLite::CrsMatrix<ComplexOrRealType>& b,
                        const int nrindex,
                        const int ncindex,
                        const int max_nnz,
                        const PsimagLite::Vector<int>::Type& rindex,
                        const PsimagLite::Vector<int>::Type& cindex,
                        PsimagLite::CrsMatrix<ComplexOrRealType>& h);

int csc_nnz(const int ncol_A,
            const PsimagLite::Vector<int>::Type& acolptr );

template<typename ComplexOrRealType>
void csc_matmul_pre( const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<ComplexOrRealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<ComplexOrRealType>& xout );

template<typename ComplexOrRealType>
void csc_matmul_post(const char trans_A, 
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<ComplexOrRealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<ComplexOrRealType>& xout);

template<typename ComplexOrRealType>
void csc_kron_mult_method(const int imethod,
                          const int nrow_A,
                          const int ncol_A,
                          const PsimagLite::Vector<int>::Type& acolptr,
                          const PsimagLite::Vector<int>::Type& arow,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                          const int nrow_B,
                          const int ncol_B,
                          const PsimagLite::Vector<int>::Type& bcolptr,
                          const PsimagLite::Vector<int>::Type& brow,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,
                          const PsimagLite::Matrix<ComplexOrRealType>& yin,
                          PsimagLite::Matrix<ComplexOrRealType>& xout );

template<typename ComplexOrRealType>
void csc_kron_mult(const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& acolptr,
                   const PsimagLite::Vector<int>::Type& arow,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& bcolptr,
                   const PsimagLite::Vector<int>::Type& brow,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,
                   const PsimagLite::Matrix<ComplexOrRealType>& yin,
                   PsimagLite::Matrix<ComplexOrRealType>& xout );

template<typename ComplexOrRealType>
void coord2csr(const int nrow_A,
               const int ncol_A,
               const int nnz,
               const int ilist[],
               const int jlist[],
               const ComplexOrRealType alist[],
               int arowptr[],
               int acol[],
               ComplexOrRealType aval[] );

template<typename ComplexOrRealType>
void den_csr_kron_mult_method(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::Matrix<ComplexOrRealType>& a_,
                              const PsimagLite::CrsMatrix<ComplexOrRealType>& b,
                              const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                              SizeType offsetY,
                              typename PsimagLite::Vector<ComplexOrRealType>::Type& xout_,
                              SizeType offsetX);

void den_copymat( const int nrow, 
                  const int ncol,
                  const int asrc_[],
                  int bdest_[] );

template<typename ComplexOrRealType>
void den_zeros( const int nrow_A,
                const int ncol_A,
                PsimagLite::Matrix<ComplexOrRealType>& a_);

template<typename ComplexOrRealType>
void den_transpose( const int nrow_A, 
                    const int ncol_A,
                    const ComplexOrRealType a_[],
                    ComplexOrRealType at_[] );

template<typename ComplexOrRealType>
void den_gen_matrix( const int nrow_A, 
                     const int ncol_A,
                     const ComplexOrRealType threshold,
                     PsimagLite::Matrix<ComplexOrRealType>& a_);

template<typename ComplexOrRealType>
void den_matmul_pre(const char trans_A,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<ComplexOrRealType>& a_,
                    const int nrow_Y,
                    const int ncol_Y,
                    const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                    const int nrow_X,
                    const int ncol_X,
                    PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout);

template<typename ComplexOrRealType>
void den_matmul_post(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<ComplexOrRealType>& a_,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const ComplexOrRealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<ComplexOrRealType>& xout);

template<typename ComplexOrRealType>
void den_kron_submatrix(
        const int nrow_A,
        const int ncol_A,
        const PsimagLite::Matrix<ComplexOrRealType>& a_,
        const int nrow_B,
        const int ncol_B,
        const PsimagLite::Matrix<ComplexOrRealType>& b_,
        const int nrindex,
        const int ncindex,
        const PsimagLite::Vector<int>::Type& rindex,
        const PsimagLite::Vector<int>::Type& cindex,
        PsimagLite::Matrix<ComplexOrRealType>& c_ );

template<typename ComplexOrRealType>
void den_kron_mult_method(const int imethod,
                          const char transA,
                          const char transB,
                          const PsimagLite::Matrix<ComplexOrRealType>& a_,
                          const PsimagLite::Matrix<ComplexOrRealType>& b_,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                          SizeType offsetY ,
                          typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
                          SizeType offsetX);

template<typename ComplexOrRealType>
int den_nnz(const PsimagLite::Matrix<ComplexOrRealType>&);

template<typename ComplexOrRealType>
bool den_is_eye(const PsimagLite::Matrix<ComplexOrRealType>&);

template<typename ComplexOrRealType>
bool den_is_zeros(const PsimagLite::Matrix<ComplexOrRealType>&);


template<typename ComplexOrRealType>
void den_kron_form( const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<ComplexOrRealType>& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<ComplexOrRealType>& b_,
                    PsimagLite::Matrix<ComplexOrRealType>& c_);


template<typename ComplexOrRealType>
void den_submatrix( const int nrow_A, 
                    const int ncol_A,
                    const PsimagLite::Matrix<ComplexOrRealType>& a_,
                    const int nrindex,
                    const int ncindex,
                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,
                    PsimagLite::Matrix<ComplexOrRealType>& c_);

template<typename ComplexOrRealType>
void den_eye(const int nrow_A,
             const int ncol_A,
             PsimagLite::Matrix<ComplexOrRealType>& c_);
#endif

