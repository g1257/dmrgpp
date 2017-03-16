#include "util.h"
int csr_nnz( const int nrow_A,  
             const PsimagLite::Vector<int>::Type& arowptr)
{
 return( arowptr[nrow_A] );
}
 
