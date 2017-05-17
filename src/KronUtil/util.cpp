#include "estimate_kron_cost.cpp"
#include "csr_den_kron_mult.cpp"
#include "csr_kron_mult.cpp"
#include "csr_nnz.cpp"
#include "csr_eye.cpp"
#include "csr_is_eye.cpp"
#include "csr_is_zeros.cpp"
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

template
void estimate_kron_cost
<double>(const int nrow_A,
                         const int ncol_A,
                         const int nnz_A,
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B,
                         double *p_kron_nnz,
                         double *p_kron_flops,
                         int *p_imethod );

template
void csr_den_kron_mult_method<double>(const int imethod,
                              const char transA,
                              const char transB,

                              const PsimagLite::CrsMatrix<double>& a_,
                              const PsimagLite::Matrix<double>& b_,

                              const typename PsimagLite::Vector<double>::Type& yin_,
                              SizeType offsetY ,
                              typename PsimagLite::Vector<double>::Type& xout_,
                              SizeType offsetX);

template
int csr_nnz<double>(const PsimagLite::CrsMatrix<double>&);

template
bool csr_is_eye<double>(const PsimagLite::CrsMatrix<double>&);

template
bool csr_is_zeros<double>(const PsimagLite::CrsMatrix<double>&);

template
void csr_transpose<double>(const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const double aval[],
                   int atrowptr[],
                   int atcol[],
                   double atval[] );

template
void csr_kron_mult_method<double>(const int imethod,
                          const char transA,
                          const char transB,

                          const PsimagLite::CrsMatrix<double>& a,
                          const PsimagLite::CrsMatrix<double>& b,

                          const PsimagLite::MatrixNonOwned<const double>& yin,
                          PsimagLite::MatrixNonOwned<double>& xout);



template
void csr_matmul_post<double>(const char trans_A,
                     const PsimagLite::CrsMatrix<double>&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);

template
void csr_matmul_pre<double>(const char trans_A,
                     const PsimagLite::CrsMatrix<double>&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);

template
void csr_submatrix<double>(const PsimagLite::CrsMatrix<double>& a,
                   const int nrow_B,
                   const int ncol_B,
                   const int max_nnz,
                   const PsimagLite::Vector<int>::Type& rindex,
                   const PsimagLite::Vector<int>::Type& cindex,
                   PsimagLite::CrsMatrix<double>& b);


template
void csr_eye<double>(const int nrow_B,
             const int ncol_B,
             PsimagLite::CrsMatrix<double>& b);

template
void csr_kron_submatrix<double>(const PsimagLite::CrsMatrix<double>& a,
                        const PsimagLite::CrsMatrix<double>& b,
                        const int nrindex,
                        const int ncindex,
                        const int max_nnz,
                        const PsimagLite::Vector<int>::Type& rindex,
                        const PsimagLite::Vector<int>::Type& cindex,
                        PsimagLite::CrsMatrix<double>& h);

template
void csc_matmul_pre<double>(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const typename PsimagLite::Vector<double>::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<double>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<double>& xout );

template
void csc_matmul_post<double>(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const typename PsimagLite::Vector<double>::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<double>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<double>& xout);

template
void csc_kron_mult_method<double>(const int imethod,
                          const int nrow_A,
                          const int ncol_A,
                          const PsimagLite::Vector<int>::Type& acolptr,
                          const PsimagLite::Vector<int>::Type& arow,
                          const typename PsimagLite::Vector<double>::Type& aval,
                          const int nrow_B,
                          const int ncol_B,
                          const PsimagLite::Vector<int>::Type& bcolptr,
                          const PsimagLite::Vector<int>::Type& brow,
                          const typename PsimagLite::Vector<double>::Type& bval,
                          const PsimagLite::Matrix<double>& yin,
                          PsimagLite::Matrix<double>& xout );

template
void csc_kron_mult<double>(const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& acolptr,
                   const PsimagLite::Vector<int>::Type& arow,
                   const typename PsimagLite::Vector<double>::Type& aval,
                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& bcolptr,
                   const PsimagLite::Vector<int>::Type& brow,
                   const typename PsimagLite::Vector<double>::Type& bval,
                   const PsimagLite::Matrix<double>& yin,
                   PsimagLite::Matrix<double>& xout );

template
void coord2csr<double>(const int nrow_A,
               const int ncol_A,
               const int nnz,
               const int ilist[],
               const int jlist[],
               const double alist[],
               int arowptr[],
               int acol[],
               double aval[] );

template
void den_csr_kron_mult_method<double>(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::Matrix<double>& a_,
                              const PsimagLite::CrsMatrix<double>& b,
                              const typename PsimagLite::Vector<double>::Type& yin,
                              SizeType offsetY,
                              typename PsimagLite::Vector<double>::Type& xout_,
                              SizeType offsetX);

template
void den_zeros<double>(const int nrow_A,
                const int ncol_A,
                PsimagLite::Matrix<double>& a_);

template
void den_transpose<double>(const int nrow_A,
                    const int ncol_A,
                    const double a_[],
                    double at_[] );

//template
//void den_gen_matrix<double>(const int nrow_A,
//                     const int ncol_A,
//                     const double threshold,
//                     PsimagLite::Matrix<double>& a_);

template
void den_matmul_pre<double>(const char trans_A,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<double>& a_,
                    const int nrow_Y,
                    const int ncol_Y,
                    const PsimagLite::MatrixNonOwned<const double>& yin,
                    const int nrow_X,
                    const int ncol_X,
                    PsimagLite::MatrixNonOwned<double>& xout);

template
void den_matmul_post<double>(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<double>& a_,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const double>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<double>& xout);

template
void den_kron_submatrix<double>(const int nrow_A,
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

template
void den_kron_mult_method<double>(const int imethod,
                          const char transA,
                          const char transB,
                          const PsimagLite::Matrix<double>& a_,
                          const PsimagLite::Matrix<double>& b_,
                          const typename PsimagLite::Vector<double>::Type& yin,
                          SizeType offsetY ,
                          typename PsimagLite::Vector<double>::Type& xout,
                          SizeType offsetX);

template
int den_nnz<double>(const PsimagLite::Matrix<double>&);

template
bool den_is_eye<double>(const PsimagLite::Matrix<double>&);

template
bool den_is_zeros<double>(const PsimagLite::Matrix<double>&);


template
void den_kron_form<double>(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<double>& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<double>& b_,
                    PsimagLite::Matrix<double>& c_);


template
void den_submatrix<double>(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<double>& a_,
                    const int nrindex,
                    const int ncindex,
                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,
                    PsimagLite::Matrix<double>& c_);

template
void den_eye<double>(const int nrow_A,
             const int ncol_A,
             PsimagLite::Matrix<double>& c_);

