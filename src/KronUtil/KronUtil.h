#ifndef KRON_UTIL_HEADER_H
#define KRON_UTIL_HEADER_H
#include <complex>
#include "Vector.h"
#include "Matrix.h"

extern
void csr_kron_mult(const char transA,
                   const char transB,
                   const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& arowptr,
                   const PsimagLite::Vector<int>::Type& acol,
                   const PsimagLite::Vector<double>::Type& aval,

                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& browptr,
                   const PsimagLite::Vector<int>::Type& bcol,
                   const PsimagLite::Vector<double>::Type& bval,

                   const double* yin,
                   double* xout);

extern
void csr_kron_mult(const char transA,
                   const char transB,
                   const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& arowptr,
                   const PsimagLite::Vector<int>::Type& acol,
                   const typename PsimagLite::Vector<std::complex<double> >::Type& aval,

                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& browptr,
                   const PsimagLite::Vector<int>::Type& bcol,
                   const typename PsimagLite::Vector<std::complex<double> >::Type& bval,

                   const std::complex<double>* yin,
                   std::complex<double>* xout);
//-----------------------------------------------------------------------------------

extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<double>& a_,

                       const int nrow_B,
                       const int ncol_B,
                       const PsimagLite::Vector<int>::Type& browptr,
                       const PsimagLite::Vector<int>::Type& bcol,
                       const PsimagLite::Vector<double>::Type& bval,

                       const double* yin_,
                       double* xout_);

extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<std::complex<double> >& a_,

                       const int nrow_B,
                       const int ncol_B,
                       const PsimagLite::Vector<int>::Type& browptr,
                       const PsimagLite::Vector<int>::Type& bcol,
                       const PsimagLite::Vector<std::complex<double> >::Type& bval,

                       const std::complex<double>* yin_,
                       std::complex<double>* xout_);

//-----------------------------------------------------------------------------------

extern
void den_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::Matrix<double>& a_,
                   const PsimagLite::Matrix<double>& b_,

                   const double* yin,
                   double* xout);

extern
void den_kron_mult(const char transA,
                   const char transB,

                   const PsimagLite::Matrix<std::complex<double> >& a_,
                   const PsimagLite::Matrix<std::complex<double> >& b_,

                   const std::complex<double>* yin,
                   std::complex<double>* xout);

//-----------------------------------------------------------------------------------

extern
void csr_den_kron_mult( const char transA,
                        const char transB,

                        const int nrow_A,
                        const int ncol_A,
                        const PsimagLite::Vector<int>::Type& arowptr,
                        const PsimagLite::Vector<int>::Type& acol,
                        const PsimagLite::Vector<double>::Type& aval,

                        const PsimagLite::Matrix<double>& b_,

                        const double* yin_,
                        double* xout_);

extern
void csr_den_kron_mult( const char transA,
                        const char transB,

                        const int nrow_A,
                        const int ncol_A,
                        const PsimagLite::Vector<int>::Type& arowptr,
                        const PsimagLite::Vector<int>::Type& acol,
                        const PsimagLite::Vector<std::complex<double> >::Type& aval,

                        const PsimagLite::Matrix<std::complex<double> >& b_,

                        const std::complex<double>* yin_,
                        std::complex<double>* xout_);

#endif

