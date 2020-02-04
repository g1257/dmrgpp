#ifndef KRON_UTIL_HEADER_H
#define KRON_UTIL_HEADER_H
#include <complex>
#include "Vector.h"
#include "Matrix.h"
#include "CrsMatrix.h"
#include "GemmR.h"

template<typename ComplexOrRealType>
void csr_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>& a,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>& b,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                   SizeType offsetY,
                   typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
                   SizeType offsetX,
                   const typename PsimagLite::Real<ComplexOrRealType>::Type);

//-----------------------------------------------------------------------------------

template<typename ComplexOrRealType>
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<ComplexOrRealType>& a_,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
	                   SizeType offsetY,
	                   typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
	                   SizeType offsetX,
                       const typename PsimagLite::Real<ComplexOrRealType>::Type,
                       PsimagLite::GemmR<ComplexOrRealType>&);

//-----------------------------------------------------------------------------------

template<typename ComplexOrRealType>
void den_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::Matrix<ComplexOrRealType>& a_,
                   const PsimagLite::Matrix<ComplexOrRealType>& b_,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                   SizeType offsetY,
                   typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
                   SizeType offsetX,
                   const typename PsimagLite::Real<ComplexOrRealType>::Type,
                   PsimagLite::GemmR<ComplexOrRealType>&);

//-----------------------------------------------------------------------------------

template<typename ComplexOrRealType>
void csr_den_kron_mult( const char transA,
                        const char transB,
                        const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                        const PsimagLite::Matrix<ComplexOrRealType>& b_,
                        const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
	                    SizeType offsetY,
	                    typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
	                    SizeType offsetX,
                        const typename PsimagLite::Real<ComplexOrRealType>::Type,
                        PsimagLite::GemmR<ComplexOrRealType>&);
#endif

