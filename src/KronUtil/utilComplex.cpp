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

template
void estimate_kron_cost
<std::complex<double> >(const int nrow_A,
                         const int ncol_A,
                         const int nnz_A,
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B,
                         std::complex<double>  *p_kron_nnz,
                         std::complex<double>  *p_kron_flops,
                         int *p_imethod );

template
void csr_den_kron_mult_method<std::complex<double> >(const int imethod,
                              const char transA,
                              const char transB,

                              const PsimagLite::CrsMatrix<std::complex<double> >& a_,
                              const PsimagLite::Matrix<std::complex<double> >& b_,

                              const PsimagLite::Vector<std::complex<double> >::Type& yin_,
                              SizeType offsetY ,
                              PsimagLite::Vector<std::complex<double> >::Type& xout_,
                              SizeType offsetX);

template
bool csr_is_eye<std::complex<double> >(const PsimagLite::CrsMatrix<std::complex<double> >&);

template
void csr_transpose<std::complex<double> >(const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const std::complex<double>  aval[],
                   int atrowptr[],
                   int atcol[],
                   std::complex<double>  atval[] );

template
void csr_kron_mult_method<std::complex<double> >(const int imethod,
                          const char transA,
                          const char transB,

                          const PsimagLite::CrsMatrix<std::complex<double> >& a,
                          const PsimagLite::CrsMatrix<std::complex<double> >& b,

                          const PsimagLite::MatrixNonOwned<const std::complex<double> >& yin,
                          PsimagLite::MatrixNonOwned<std::complex<double> >& xout);



template
void csr_matmul_post<std::complex<double> >(const char trans_A,
                     const PsimagLite::CrsMatrix<std::complex<double> >&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const std::complex<double> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<std::complex<double> >& xout);

template
void csr_matmul_pre<std::complex<double> >(const char trans_A,
                     const PsimagLite::CrsMatrix<std::complex<double> >&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const std::complex<double> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<std::complex<double> >& xout);

template
void csr_submatrix<std::complex<double> >(const PsimagLite::CrsMatrix<std::complex<double> >& a,
                   const int nrow_B,
                   const int ncol_B,
                   const int max_nnz,
                   const PsimagLite::Vector<int>::Type& rindex,
                   const PsimagLite::Vector<int>::Type& cindex,
                   PsimagLite::CrsMatrix<std::complex<double> >& b);


template
void csr_eye<std::complex<double> >(const int nrow_B,
             const int ncol_B,
             PsimagLite::CrsMatrix<std::complex<double> >& b);

template
void csr_kron_submatrix<std::complex<double> >(const PsimagLite::CrsMatrix<std::complex<double> >& a,
                        const PsimagLite::CrsMatrix<std::complex<double> >& b,
                        const int nrindex,
                        const int ncindex,
                        const int max_nnz,
                        const PsimagLite::Vector<int>::Type& rindex,
                        const PsimagLite::Vector<int>::Type& cindex,
                        PsimagLite::CrsMatrix<std::complex<double> >& h);

template
void csc_matmul_pre<std::complex<double> >(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<std::complex<double> >::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<std::complex<double> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<std::complex<double> >& xout );

template
void csc_matmul_post<std::complex<double> >(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<std::complex<double> >::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<std::complex<double> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<std::complex<double> >& xout);

template
void csc_kron_mult_method<std::complex<double> >(const int imethod,
                          const int nrow_A,
                          const int ncol_A,
                          const PsimagLite::Vector<int>::Type& acolptr,
                          const PsimagLite::Vector<int>::Type& arow,
                          const PsimagLite::Vector<std::complex<double> >::Type& aval,
                          const int nrow_B,
                          const int ncol_B,
                          const PsimagLite::Vector<int>::Type& bcolptr,
                          const PsimagLite::Vector<int>::Type& brow,
                          const PsimagLite::Vector<std::complex<double> >::Type& bval,
                          const PsimagLite::Matrix<std::complex<double> >& yin,
                          PsimagLite::Matrix<std::complex<double> >& xout );

template
void csc_kron_mult<std::complex<double> >(const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& acolptr,
                   const PsimagLite::Vector<int>::Type& arow,
                   const PsimagLite::Vector<std::complex<double> >::Type& aval,
                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& bcolptr,
                   const PsimagLite::Vector<int>::Type& brow,
                   const PsimagLite::Vector<std::complex<double> >::Type& bval,
                   const PsimagLite::Matrix<std::complex<double> >& yin,
                   PsimagLite::Matrix<std::complex<double> >& xout );

template
void coord2csr<std::complex<double> >(const int nrow_A,
               const int ncol_A,
               const int nnz,
               const int ilist[],
               const int jlist[],
               const std::complex<double>  alist[],
               int arowptr[],
               int acol[],
               std::complex<double>  aval[] );

template
void den_csr_kron_mult_method<std::complex<double> >(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::Matrix<std::complex<double> >& a_,
                              const PsimagLite::CrsMatrix<std::complex<double> >& b,
                              const PsimagLite::Vector<std::complex<double> >::Type& yin,
                              SizeType offsetY,
                              PsimagLite::Vector<std::complex<double> >::Type& xout_,
                              SizeType offsetX);

template
void den_zeros<std::complex<double> >(const int nrow_A,
                const int ncol_A,
                PsimagLite::Matrix<std::complex<double> >& a_);

template
void den_transpose<std::complex<double> >(const int nrow_A,
                    const int ncol_A,
                    const std::complex<double>  a_[],
                    std::complex<double>  at_[] );

template
void den_gen_matrix<std::complex<double> >(const int nrow_A,
                     const int ncol_A,
                     const double& threshold,
                     PsimagLite::Matrix<std::complex<double> >& a_);

template
void den_matmul_pre<std::complex<double> >(const char trans_A,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<std::complex<double> >& a_,
                    const int nrow_Y,
                    const int ncol_Y,
                    const PsimagLite::MatrixNonOwned<const std::complex<double> >& yin,
                    const int nrow_X,
                    const int ncol_X,
                    PsimagLite::MatrixNonOwned<std::complex<double> >& xout);

template
void den_matmul_post<std::complex<double> >(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<std::complex<double> >& a_,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const std::complex<double> >& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<std::complex<double> >& xout);

template
void den_kron_submatrix<std::complex<double> >(const int nrow_A,
        const int ncol_A,
        const PsimagLite::Matrix<std::complex<double> >& a_,
        const int nrow_B,
        const int ncol_B,
        const PsimagLite::Matrix<std::complex<double> >& b_,
        const int nrindex,
        const int ncindex,
        const PsimagLite::Vector<int>::Type& rindex,
        const PsimagLite::Vector<int>::Type& cindex,
        PsimagLite::Matrix<std::complex<double> >& c_ );

template
void den_kron_mult_method<std::complex<double> >(const int imethod,
                          const char transA,
                          const char transB,
                          const PsimagLite::Matrix<std::complex<double> >& a_,
                          const PsimagLite::Matrix<std::complex<double> >& b_,
                          const PsimagLite::Vector<std::complex<double> >::Type& yin,
                          SizeType offsetY ,
                          PsimagLite::Vector<std::complex<double> >::Type& xout,
                          SizeType offsetX);

template
int den_nnz<std::complex<double> >(const PsimagLite::Matrix<std::complex<double> >&);

template
bool den_is_eye<std::complex<double> >(const PsimagLite::Matrix<std::complex<double> >&);

template
bool den_is_zeros<std::complex<double> >(const PsimagLite::Matrix<std::complex<double> >&);


template
void den_kron_form<std::complex<double> >(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<std::complex<double> >& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<std::complex<double> >& b_,
                    PsimagLite::Matrix<std::complex<double> >& c_);


template
void den_submatrix<std::complex<double> >(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<std::complex<double> >& a_,
                    const int nrindex,
                    const int ncindex,
                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,
                    PsimagLite::Matrix<std::complex<double> >& c_);

template
void den_eye<std::complex<double> >(const int nrow_A,
             const int ncol_A,
             PsimagLite::Matrix<std::complex<double> >& c_);

