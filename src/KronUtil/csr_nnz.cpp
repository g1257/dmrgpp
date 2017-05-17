#include "util.h"

template<typename ComplexOrRealType>
int csr_nnz(const PsimagLite::CrsMatrix<ComplexOrRealType>& a)
{
	return a.getRowPtr(a.row());
}
 
