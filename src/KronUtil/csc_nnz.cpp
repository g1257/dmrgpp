#include "util.h"
int csc_nnz( const int ncol_A, 
             const PsimagLite::Vector<int>::Type& acolptr)
{
/*
 * ---------------------------------
 * return the number of nonzeros
 * matrix in compress COLUMN format
 * ---------------------------------
 */
 return( acolptr[ncol_A] );
}
