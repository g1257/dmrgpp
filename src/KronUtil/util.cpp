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
#include "den_kron_form_general.cpp"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

template
void estimate_kron_cost
<RealType>(const int nrow_A,
                         const int ncol_A,
                         const int nnz_A,
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B,
                         RealType *p_kron_nnz,
                         RealType *p_kron_flops,
                         int *p_imethod,
                         const RealType dense_flop_discount);

template
void csr_den_kron_mult_method<RealType>(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::CrsMatrix<RealType>& a_,
                              const PsimagLite::Matrix<RealType>& b_,
                              const PsimagLite::Vector<RealType>::Type& yin_,
                              SizeType offsetY ,
                              PsimagLite::Vector<RealType>::Type& xout_,
                              SizeType offsetX,
                              PsimagLite::GemmR<RealType>&);


template
bool csr_is_eye<RealType>(const PsimagLite::CrsMatrix<RealType>&);

template
void csr_transpose<RealType>(const int nrow_A,
                   const int ncol_A,
                   const int arowptr[],
                   const int acol[],
                   const RealType aval[],
                   int atrowptr[],
                   int atcol[],
                   RealType atval[] );

template
void csr_kron_mult_method<RealType>(const int imethod,
                          const char transA,
                          const char transB,

                          const PsimagLite::CrsMatrix<RealType>& a,
                          const PsimagLite::CrsMatrix<RealType>& b,

                          const PsimagLite::MatrixNonOwned<const RealType>& yin,
                          PsimagLite::MatrixNonOwned<RealType>& xout);



template
void csr_matmul_post<RealType>(const char trans_A,
                     const PsimagLite::CrsMatrix<RealType>&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const RealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<RealType>& xout);

template
void csr_matmul_pre<RealType>(const char trans_A,
                     const PsimagLite::CrsMatrix<RealType>&,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const RealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<RealType>& xout);

template
void csr_submatrix<RealType>(const PsimagLite::CrsMatrix<RealType>& a,
                   const int nrow_B,
                   const int ncol_B,
                   const int max_nnz,
                   const PsimagLite::Vector<int>::Type& rindex,
                   const PsimagLite::Vector<int>::Type& cindex,
                   PsimagLite::CrsMatrix<RealType>& b);


template
void csr_eye<RealType>(const int nrow_B,
             const int ncol_B,
             PsimagLite::CrsMatrix<RealType>& b);

template
void csr_kron_submatrix<RealType>(const PsimagLite::CrsMatrix<RealType>& a,
                        const PsimagLite::CrsMatrix<RealType>& b,
                        const int nrindex,
                        const int ncindex,
                        const int max_nnz,
                        const PsimagLite::Vector<int>::Type& rindex,
                        const PsimagLite::Vector<int>::Type& cindex,
                        PsimagLite::CrsMatrix<RealType>& h);

template
void csc_matmul_pre<RealType>(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<RealType>::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<RealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<RealType>& xout );

template
void csc_matmul_post<RealType>(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Vector<int>::Type& acolptr,
                     const PsimagLite::Vector<int>::Type& arow,
                     const PsimagLite::Vector<RealType>::Type& aval,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::Matrix<RealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::Matrix<RealType>& xout);

template
void csc_kron_mult_method<RealType>(const int imethod,
                          const int nrow_A,
                          const int ncol_A,
                          const PsimagLite::Vector<int>::Type& acolptr,
                          const PsimagLite::Vector<int>::Type& arow,
                          const PsimagLite::Vector<RealType>::Type& aval,
                          const int nrow_B,
                          const int ncol_B,
                          const PsimagLite::Vector<int>::Type& bcolptr,
                          const PsimagLite::Vector<int>::Type& brow,
                          const PsimagLite::Vector<RealType>::Type& bval,
                          const PsimagLite::Matrix<RealType>& yin,
                          PsimagLite::Matrix<RealType>& xout );

template
void csc_kron_mult<RealType>(const int nrow_A,
                             const int ncol_A,
                             const PsimagLite::Vector<int>::Type& acolptr,
                             const PsimagLite::Vector<int>::Type& arow,
                             const PsimagLite::Vector<RealType>::Type& aval,
                             const int nrow_B,
                             const int ncol_B,
                             const PsimagLite::Vector<int>::Type& bcolptr,
                             const PsimagLite::Vector<int>::Type& brow,
                             const PsimagLite::Vector<RealType>::Type& bval,
                             const PsimagLite::Matrix<RealType>& yin,
                             PsimagLite::Matrix<RealType>& xout,
                             const RealType);

template
void coord2csr<RealType>(const int nrow_A,
               const int ncol_A,
               const int nnz,
               const int ilist[],
               const int jlist[],
               const RealType alist[],
               int arowptr[],
               int acol[],
               RealType aval[] );

template
void den_csr_kron_mult_method<RealType>(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::Matrix<RealType>& a_,
                              const PsimagLite::CrsMatrix<RealType>& b,
                              const PsimagLite::Vector<RealType>::Type& yin,
                              SizeType offsetY,
                              PsimagLite::Vector<RealType>::Type& xout_,
                              SizeType offsetX,
                              PsimagLite::GemmR<RealType>&);

template
void den_zeros<RealType>(const int nrow_A,
                const int ncol_A,
                PsimagLite::Matrix<RealType>& a_);

template
void den_transpose<RealType>(const int nrow_A,
                    const int ncol_A,
                    const RealType a_[],
                    RealType at_[] );

template
void den_gen_matrix<RealType>(const int nrow_A,
                     const int ncol_A,
                     const RealType& threshold,
                     PsimagLite::Matrix<RealType>& a_);

template
void den_matmul_pre<RealType>(const char trans_A,
                    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<RealType>& a_,
                    const int nrow_Y,
                    const int ncol_Y,
                    const PsimagLite::MatrixNonOwned<const RealType>& yin,
                    const int nrow_X,
                    const int ncol_X,
                    PsimagLite::MatrixNonOwned<RealType>& xout,
                    PsimagLite::GemmR<RealType>&);

template
void den_matmul_post<RealType>(const char trans_A,
                     const int nrow_A,
                     const int ncol_A,
                     const PsimagLite::Matrix<RealType>& a_,
                     const int nrow_Y,
                     const int ncol_Y,
                     const PsimagLite::MatrixNonOwned<const RealType>& yin,
                     const int nrow_X,
                     const int ncol_X,
                     PsimagLite::MatrixNonOwned<RealType>& xout,
                     PsimagLite::GemmR<RealType>&);

template
void den_kron_submatrix<RealType>(const int nrow_A,
        const int ncol_A,
        const PsimagLite::Matrix<RealType>& a_,
        const int nrow_B,
        const int ncol_B,
        const PsimagLite::Matrix<RealType>& b_,
        const int nrindex,
        const int ncindex,
        const PsimagLite::Vector<int>::Type& rindex,
        const PsimagLite::Vector<int>::Type& cindex,
        PsimagLite::Matrix<RealType>& c_ );

template
void den_kron_mult_method<RealType>(const int imethod,
                          const char transA,
                          const char transB,
                          const PsimagLite::Matrix<RealType>& a_,
                          const PsimagLite::Matrix<RealType>& b_,
                          const PsimagLite::Vector<RealType>::Type& yin,
                          SizeType offsetY ,
                          PsimagLite::Vector<RealType>::Type& xout,
                          SizeType offsetX,
                          PsimagLite::GemmR<RealType>&);

template
int den_nnz<RealType>(const PsimagLite::Matrix<RealType>&);

template
bool den_is_eye<RealType>(const PsimagLite::Matrix<RealType>&);

template
bool den_is_zeros<RealType>(const PsimagLite::Matrix<RealType>&);


template
void den_kron_form<RealType>(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<RealType>& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<RealType>& b_,
                    PsimagLite::Matrix<RealType>& c_);


template
void den_kron_form_general<RealType>(
		    const char transA,
		    const char transB,
		    const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<RealType>& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<RealType>& b_,
                    PsimagLite::Matrix<RealType>& c_);

template
void den_submatrix<RealType>(const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Matrix<RealType>& a_,
                    const int nrindex,
                    const int ncindex,
                    const PsimagLite::Vector<int>::Type& rindex,
                    const PsimagLite::Vector<int>::Type& cindex,
                    PsimagLite::Matrix<RealType>& c_);

template
void den_eye<RealType>(const int nrow_A,
             const int ncol_A,
             PsimagLite::Matrix<RealType>& c_);

