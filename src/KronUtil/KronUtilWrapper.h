#ifndef KRON_UTIL_WRAPPER_H
#define KRON_UTIL_WRAPPER_H

#ifndef DO_NOT_USE_KRON_UTIL
#include "KronUtil.h"
#else
#include "Matrix.h"
#include "ProgramGlobals.h"

template <typename ComplexOrRealType>
void csr_kron_mult(const char transA,
    const char transB,
    const PsimagLite::CrsMatrix<ComplexOrRealType>&,
    const PsimagLite::CrsMatrix<ComplexOrRealType>&,
    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
    SizeType offsetY,
    typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
    SizeType offsetX)

{
	PsimagLite::String msg("csr_kron_mult: please #undefine DO_NOT_USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template <typename ComplexOrRealType>
void csr_den_kron_mult(const char transA,
    const char transB,
    const PsimagLite::CrsMatrix<ComplexOrRealType>&,
    const PsimagLite::Matrix<ComplexOrRealType>& b_,
    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
    SizeType offsetY,
    typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
    SizeType offsetX)
{
	PsimagLite::String msg("csr_den_kron_mult: please #undefine DO_NOT_USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template <typename ComplexOrRealType>
void den_csr_kron_mult(const char transA,
    const char transB,
    const PsimagLite::Matrix<ComplexOrRealType>& a_,
    const PsimagLite::CrsMatrix<ComplexOrRealType>&,
    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
    SizeType offsetY,
    typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
    SizeType offsetX)
{
	PsimagLite::String msg("den_csr_kron_mult: please #undefine DO_NOT_USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template <typename ComplexOrRealType>
void den_kron_mult(const char transA,
    const char transB,
    const PsimagLite::Matrix<ComplexOrRealType>& a_,
    const PsimagLite::Matrix<ComplexOrRealType>& b_,
    const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
    SizeType offsetY,
    typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
    SizeType offsetX)
{
	PsimagLite::String msg("den_kron_mult: please #undefine DO_NOT_USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

#endif

#endif // KRON_UTIL_WRAPPER_H
