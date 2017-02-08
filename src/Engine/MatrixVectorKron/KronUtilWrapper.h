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
                   const int arowptr[],
                   const int acol[],
                   const ComplexOrRealType aval[],

                   const int nrow_B,
                   const int ncol_B,
                   const int browptr[],
                   const int bcol[],
                   const ComplexOrRealType bval[],

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
                       const int arowptr[],
                       const int acol[],
                       const ComplexOrRealType aval[],

                       const int nrow_B,
                       const int ncol_B,
                       const ComplexOrRealType b_[],

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
                       const ComplexOrRealType a_[],

                       const int nrow_B,
                       const int ncol_B,
                       const int browptr[],
                       const int bcol[],
                       const ComplexOrRealType bval[],

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
                   const ComplexOrRealType a_[],

                   const int nrow_B,
                   const int ncol_B,
                   const ComplexOrRealType b_[],

                   const ComplexOrRealType yin[],
                   ComplexOrRealType xout[])
{
	PsimagLite::String msg("den_kron_mult: please #define USE_KRON_UTIL");
	msg += " and link against libkronutil\n";
	throw PsimagLite::RuntimeError(msg);
}

#endif

#endif // KRON_UTIL_WRAPPER_H
