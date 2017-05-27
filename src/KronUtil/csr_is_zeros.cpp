#include "util.h"

template<typename ComplexOrRealType>
bool csr_is_zeros(const PsimagLite::CrsMatrix<ComplexOrRealType>& a)
{
	// ----------------------------------------------------
	// check whether a sparse matrix is the zero matrix
	// ----------------------------------------------------
	return isZero(a, 1e-12);
}

