#include "estimate_kron_cost.cpp"
#include "csr_den_kron_mult.cpp"
#include "csr_kron_mult.cpp"
#include "csr_eye.cpp"
#include "csr_is_eye.cpp"
#include "csr_transpose.cpp"
#include "csr_matmul_post.cpp"
#include "csr_matmul_pre.cpp"
#include "csr_submatrix.cpp"
#include "csr_kron_submatrix.cpp"
#include "csc_matmul_pre.cpp"
#include "csc_matmul_post.cpp"
#include "csc_kron_mult.cpp"
#include "coord2csr.cpp"
#include "den_csr_kron_mult.cpp"
#include "den_zeros.cpp"
#include "den_transpose.cpp"
#include "den_gen_matrix.cpp"
#include "den_matmul_pre.cpp"
#include "den_matmul_post.cpp"
#include "den_kron_submatrix.cpp"
#include "den_submatrix.cpp"
#include "den_kron_mult.cpp"
#include "den_nnz.cpp"
#include "den_eye.cpp"
#include "den_is_eye.cpp"
#include "den_is_zeros.cpp"
#include "den_kron_form.cpp"
#include <complex>

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

template
void estimate_kron_cost
<std::complex<RealType> >(const int nrow_A,
                         const int ncol_A,
                         const int nnz_A,
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B,
                         std::complex<RealType>  *p_kron_nnz,
                         std::complex<RealType>  *p_kron_flops,
                         int *p_imethod );

template
void csr_den_kron_mult_method<std::complex<RealType> >(const int imethod,
                              const char transA,
                              const char transB,

                              const PsimagLite::CrsMatrix<std::complex<RealType> >& a_,
                              const PsimagLite::Matrix<std::complex<RealType> >& b_,

                              const PsimagLite::Vector<std::complex<RealType> >::Type& yin_,
                              SizeType offsetY ,
                              PsimagLite::Vector<std::complex<RealType> >::Type& xout_,
                              SizeType offsetX);

template
bool csr_is_eye<std::complex<RealType> >(const PsimagLite::CrsMatrix<std::complex<RealType> >&);

template
void csr_transpose<std::complex<RealType> >(const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const std::complex<RealType>  aval[],
                   int atrowptr[],
                   int atcol[],
                   std::complex<RealType>  atval[] );

template
void csr_kron_mult_method<std::complex<RealType> >(const int imethod,
                          const char transA,
                          const char transB,

                          const PsimagLite::CrsMatrix<std::complex<RealType> >& a,
                          const PsimagLite::CrsMatrix<std::complex<RealType> >& b,

                          const PsimagLite::MatrixNonOwned<const std::complex<RealType> >& yin,
                          PsimagLite::MatrixNonOwned<std::complex<RealType> >& xout);



template
void csr_matmul_post<std::complex<RealType> >(const char trans_A,
                     const PsimagLite::CrsMatrix<std::complex<RealType> >&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const std::complex<RealType> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<std::complex<RealType> >& xout);

template
void csr_matmul_pre<std::complex<RealType> >(const char trans_A,
                     const PsimagLite::CrsMatrix<std::complex<RealType> >&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const std::complex<RealType> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<std::complex<RealType> >& xout);

template
void csr_submatrix<std::complex<RealType> >(const PsimagLite::CrsMatrix<std::complex<RealType> >& a,
                   const int nrow_B,
                   const int ncol_B,
                   const int max_nnz,
                   const PsimagLite::Vector<int>::Type& rindex,
                   const PsimagLite::Vector<int>::Type& cindex,
                   PsimagLite::CrsMatrix<std::complex<RealType> >& b);


template
void csr_eye<std::complex<RealType> >(const int nrow_B,
             const int ncol_B,
             PsimagLite::CrsMatrix<std::complex<RealType> >& b);

template
void csr_kron_submatrix<std::complex<RealType> >(const PsimagLite::CrsMatrix<std::complex<RealType> >& a,
                        const PsimagLite::CrsMatrix<std::complex<RealType> >& b,
                        const int nrindex,
                        const int ncindex,
                        const int max_nnz,
                        const PsimagLite::Vector<int>::Type& rindex,
                        const PsimagLite::Vector<int>::Type& cindex,
                        PsimagLite::CrsMatrix<std::complex<RealType> >& h);

template
void csc_matmul_pre<std::complex<RealType> >(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<std::complex<RealType> >::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<std::complex<RealType> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<std::complex<RealType> >& xout );

template
void csc_matmul_post<std::complex<RealType> >(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<std::complex<RealType> >::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<std::complex<RealType> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<std::complex<RealType> >& xout);

template
void csc_kron_mult_method<std::complex<RealType> >(const int imethod,
                          const int nrow_A,
                          const int ncol_A,
                          const PsimagLite::Vector<int>::Type& acolptr,
                          const PsimagLite::Vector<int>::Type& arow,
                          const PsimagLite::Vector<std::complex<RealType> >::Type& aval,
                          const int nrow_B,
                          const int ncol_B,
                          const PsimagLite::Vector<int>::Type& bcolptr,
                          const PsimagLite::Vector<int>::Type& brow,
                          const PsimagLite::Vector<std::complex<RealType> >::Type& bval,
                          const PsimagLite::Matrix<std::complex<RealType> >& yin,
                          PsimagLite::Matrix<std::complex<RealType> >& xout );

template
void csc_kron_mult<std::complex<RealType> >(const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& acolptr,
                   const PsimagLite::Vector<int>::Type& arow,
                   const PsimagLite::Vector<std::complex<RealType> >::Type& aval,
                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& bcolptr,
                   const PsimagLite::Vector<int>::Type& brow,
                   const PsimagLite::Vector<std::complex<RealType> >::Type& bval,
                   const PsimagLite::Matrix<std::complex<RealType> >& yin,
                   PsimagLite::Matrix<std::complex<RealType> >& xout );

template
void coord2csr<std::complex<RealType> >(const int nrow_A,
               const int ncol_A,
               const int nnz,
               const int ilist[],
               const int jlist[],
               const std::complex<RealType>  alist[],
               int arowptr[],
               int acol[],
               std::complex<RealType>  aval[] );

template
void den_csr_kron_mult_method<std::complex<RealType> >(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::Matrix<std::complex<RealType> >& a_,
                              const PsimagLite::CrsMatrix<std::complex<RealType> >& b,
                              const PsimagLite::Vector<std::complex<RealType> >::Type& yin,
                              SizeType offsetY,
                              PsimagLite::Vector<std::complex<RealType> >::Type& xout_,
                              SizeType offsetX);

template
void den_zeros<std::complex<RealType> >(const int nrow_A,
                const int ncol_A,
                PsimagLite::Matrix<std::complex<RealType> >& a_);

template
void den_transpose<std::complex<RealType> >(const int nrow_A,
                    const int ncol_A,
                    const std::complex<RealType>  a_[],
                    std::complex<RealType>  at_[] );

template
void den_gen_matrix<std::complex<RealType> >(const int nrow_A,
                     const int ncol_A,
                     const RealType& threshold,
                     PsimagLite::Matrix<std::complex<RealType> >& a_);

template
void den_matmul_pre<std::complex<RealType> >(const char trans_A,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<std::complex<RealType> >& a_,
                    const int nrow_Y,
                    const int ncol_Y,
                    const PsimagLite::MatrixNonOwned<const std::complex<RealType> >& yin,
                    const int nrow_X,
                    const int ncol_X,
                    PsimagLite::MatrixNonOwned<std::complex<RealType> >& xout);

template
void den_matmul_post<std::complex<RealType> >(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<std::complex<RealType> >& a_,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const std::complex<RealType> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<std::complex<RealType> >& xout);

template
void den_kron_submatrix<std::complex<RealType> >(const int nrow_A,
        const int ncol_A,
        const PsimagLite::Matrix<std::complex<RealType> >& a_,
        const int nrow_B,
        const int ncol_B,
        const PsimagLite::Matrix<std::complex<RealType> >& b_,
        const int nrindex,
        const int ncindex,
        const PsimagLite::Vector<int>::Type& rindex,
        const PsimagLite::Vector<int>::Type& cindex,
        PsimagLite::Matrix<std::complex<RealType> >& c_ );

template
void den_kron_mult_method<std::complex<RealType> >(const int imethod,
                          const char transA,
                          const char transB,
                          const PsimagLite::Matrix<std::complex<RealType> >& a_,
                          const PsimagLite::Matrix<std::complex<RealType> >& b_,
                          const PsimagLite::Vector<std::complex<RealType> >::Type& yin,
                          SizeType offsetY ,
                          PsimagLite::Vector<std::complex<RealType> >::Type& xout,
                          SizeType offsetX);

template
int den_nnz<std::complex<RealType> >(const PsimagLite::Matrix<std::complex<RealType> >&);

template
bool den_is_eye<std::complex<RealType> >(const PsimagLite::Matrix<std::complex<RealType> >&);

template
bool den_is_zeros<std::complex<RealType> >(const PsimagLite::Matrix<std::complex<RealType> >&);


template
void den_kron_form<std::complex<RealType> >(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<std::complex<RealType> >& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<std::complex<RealType> >& b_,
                    PsimagLite::Matrix<std::complex<RealType> >& c_);


template
void den_submatrix<std::complex<RealType> >(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<std::complex<RealType> >& a_,
                    const int nrindex,
                    const int ncindex,
                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,
                    PsimagLite::Matrix<std::complex<RealType> >& c_);

template
void den_eye<std::complex<RealType> >(const int nrow_A,
             const int ncol_A,
             PsimagLite::Matrix<std::complex<RealType> >& c_);

