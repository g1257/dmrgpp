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
                   const PsimagLite::Vector<double>::Type& yin,
                   SizeType offsetY,
                   PsimagLite::Vector<double>::Type& xout,
                   SizeType offsetX);

extern
void csr_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::CrsMatrix<std::complex<double> >&,
                   const PsimagLite::CrsMatrix<std::complex<double> >&,
                   const PsimagLite::Vector<std::complex<double> >::Type& yin,
                   SizeType offsetY,
                   PsimagLite::Vector<std::complex<double> >::Type& xout,
                   SizeType offsetX);

//-----------------------------------------------------------------------------------

extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<double>& a_,
                       const PsimagLite::CrsMatrix<double>&,
                       const PsimagLite::Vector<double>::Type& yin,
	                   SizeType offsetY,
	                   PsimagLite::Vector<double>::Type& xout,
	                   SizeType offsetX);
extern
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<std::complex<double> >& a_,
                       const PsimagLite::CrsMatrix<std::complex<double> >&,
                       const PsimagLite::Vector<std::complex<double> >::Type& yin,
                       SizeType offsetY,
                       PsimagLite::Vector<std::complex<double> >::Type& xout,
                       SizeType offsetX);


//-----------------------------------------------------------------------------------

extern
void den_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::Matrix<double>& a_,
                   const PsimagLite::Matrix<double>& b_,
                   const PsimagLite::Vector<double>::Type& yin,
                   SizeType offsetY,
                   PsimagLite::Vector<double>::Type& xout,
                   SizeType offsetX);

extern
void den_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::Matrix<std::complex<double> >& a_,
                   const PsimagLite::Matrix<std::complex<double> >& b_,
                   const PsimagLite::Vector<std::complex<double> >::Type& yin,
                   SizeType offsetY,
                   PsimagLite::Vector<std::complex<double> >::Type& xout,
                   SizeType offsetX);


//-----------------------------------------------------------------------------------

extern
void csr_den_kron_mult( const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<double>&,
                        const PsimagLite::Matrix<double>& b_,
                        const PsimagLite::Vector<double>::Type& yin,
	                    SizeType offsetY,
	                    PsimagLite::Vector<double>::Type& xout,
	                    SizeType offsetX);

extern
void csr_den_kron_mult( const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<std::complex<double> >&,
                        const PsimagLite::Matrix<std::complex<double> >& b_,
                        const PsimagLite::Vector<std::complex<double> >::Type& yin,
	                    SizeType offsetY,
	                    PsimagLite::Vector<std::complex<double> >::Type& xout,
	                    SizeType offsetX);

#endif

