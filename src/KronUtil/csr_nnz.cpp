#include "util.h"
int csr_nnz(const PsimagLite::CrsMatrix<double>& a)
{
	return a.getRowPtr(a.row());
}
 
