#ifndef KRON_UTIL_WRAPPER_H
#define KRON_UTIL_WRAPPER_H

#ifdef USE_KRON_UTIL
#include "../../KronUtil/KronUtil.h"
#else
#include "ProgramGlobals.h"

template<typename ComplexOrRealType>
void csr_kron_mult(const char transA,
                   const char transB,
                   const int nrow_A,
                   const int ncol_A,
                   const PsimagLite::Vector<int>::Type& arowptr,
                   const PsimagLite::Vector<int>::Type& acol,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,

                   const int nrow_B,
                   const int ncol_B,
                   const PsimagLite::Vector<int>::Type& browptr,
                   const PsimagLite::Vector<int>::Type& bcol,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,

                   const ComplexOrRealType yin[],
                   ComplexOrRealType xout[])
{
	PsimagLite::String msg("csr_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template<typename ComplexOrRealType>
void csr_den_kron_mult(const char transA,
                       const char transB,

                       const int nrow_A,
                       const int ncol_A,
                       const PsimagLite::Vector<int>::Type& arowptr,
                       const PsimagLite::Vector<int>::Type& acol,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,

                       const int nrow_B,
                       const int ncol_B,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& b_,

                       const ComplexOrRealType yin[],
                       ComplexOrRealType xout[] )
{
	PsimagLite::String msg("csr_den_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template<typename ComplexOrRealType>
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const int nrow_A,
                       const int ncol_A,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& a_,

                       const int nrow_B,
                       const int ncol_B,
                       const PsimagLite::Vector<int>::Type& browptr,
                       const PsimagLite::Vector<int>::Type& bcol,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,

                       const ComplexOrRealType yin[],
                       ComplexOrRealType xout[])
{
	PsimagLite::String msg("den_csr_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

template<typename ComplexOrRealType>
void den_kron_mult(const char transA,
                   const char transB,

                   const int nrow_A,
                   const int ncol_A,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& a_,

                   const int nrow_B,
                   const int ncol_B,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& b_,

                   const ComplexOrRealType yin[],
                   ComplexOrRealType xout[])
{
	PsimagLite::String msg("den_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

#endif

#endif // KRON_UTIL_WRAPPER_H
