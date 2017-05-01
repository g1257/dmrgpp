#include "util.h"
bool csr_is_zeros(const PsimagLite::CrsMatrix<double>& a)
{
   // ----------------------------------------------------
   // check whether a sparse matrix is the zero matrix
   // ----------------------------------------------------
        
   const double zero = 0;
   const int nrow_A = a.row();
   const int ncol_A = a.col();



   for(int ia=0; ia < nrow_A; ia++) {
     int istart = a.getRowPtr(ia);
     int iend = a.getRowPtr(ia + 1);
     for(int k=istart; k < iend; k++) {
         int ja = a.getCol(k);
 
         bool isok = ((0 <= ja) && (ja < ncol_A)); 
         assert( isok );

         double aij = a.getValue(k);
         if (isok && (aij != zero)) { return( false ); };
       };
    };

   return( true );
}
 
