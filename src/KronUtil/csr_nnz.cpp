#include "util.h"
int csr_nnz( const int nrow_A,  
             const int arowptr[] )
{
 return( arowptr[nrow_A] );
}
 
