#ifndef KRON_UTIL_HEADER_H
#define KRON_UTIL_HEADER_H
#include <complex>
#include "Vector.h"
#include "Matrix.h"
#include "CrsMatrix.h"

extern
void csr_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::CrsMatrix<double>& a,
                   const PsimagLite::CrsMatrix<double>& b,
                   const double* yin,
                   double* xout);

extern
void csr_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::CrsMatrix<std::complex<double> >&,
                   const PsimagLite::CrsMatrix<std::complex<double> >&,
                   const std::complex<double>* yin,
                   std::complex<double>* xout);
//-----------------------------------------------------------------------------------

extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<double>& a_,
                       const PsimagLite::CrsMatrix<double>&,
                       const double* yin_,
                       double* xout_);

extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<std::complex<double> >& a_,
                       const PsimagLite::CrsMatrix<std::complex<double> >&,
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
                        const PsimagLite::CrsMatrix<double>&,
                        const PsimagLite::Matrix<double>& b_,
                        const double* yin_,
                        double* xout_);

extern
void csr_den_kron_mult( const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<std::complex<double> >&,
                        const PsimagLite::Matrix<std::complex<double> >& b_,
                        const std::complex<double>* yin_,
                        std::complex<double>* xout_);

#endif

