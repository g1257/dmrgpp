#ifndef KRON_UTIL_WRAPPER_H
#define KRON_UTIL_WRAPPER_H

#ifdef USE_KRON_UTIL
#include "../../KronUtil/KronUtil.h"
#else
#include "ProgramGlobals.h"
#include "Matrix.h"

template<typename ComplexOrRealType>
void csr_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                   const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                   const PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                   SizeType offsetY,
                   PsimagLite::Vector<ComplexOrRealType>::Type& xout,
                   SizeType offsetX)

{
	PsimagLite::String msg("csr_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template<typename ComplexOrRealType>
void csr_den_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                       const PsimagLite::Matrix<ComplexOrRealType>& b_,
                       const PsimagLite::Vector<ComplexOrRealType>::Type& yin,
	                   SizeType offsetY,
	                   PsimagLite::Vector<ComplexOrRealType>::Type& xout,
	                   SizeType offsetX)
{
	PsimagLite::String msg("csr_den_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template<typename ComplexOrRealType>
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<ComplexOrRealType>& a_,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>&,
                       const PsimagLite::Vector<ComplexOrRealType>::Type& yin,
	                   SizeType offsetY,
	                   PsimagLite::Vector<ComplexOrRealType>::Type& xout,
	                   SizeType offsetX)
{
	PsimagLite::String msg("den_csr_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template<typename ComplexOrRealType>
void den_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::Matrix<ComplexOrRealType>& a_,
                   const PsimagLite::Matrix<ComplexOrRealType>& b_,
                   const PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                   SizeType offsetY,
                   PsimagLite::Vector<ComplexOrRealType>::Type& xout,
                   SizeType offsetX)
{
	PsimagLite::String msg("den_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

#endif

#endif // KRON_UTIL_WRAPPER_H
