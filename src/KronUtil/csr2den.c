#include "util.h"
void csr2den(  const int nrow_A,
               const int ncol_A, 
               const int arowptr[],
               const int acol[],
               const double aval[],

               double a_[] )
{
#define A(ia,ja) a_[(ia) + (ja)*nrow_A]
/*
 * ----------------------------------------------------------
 * convert from compress sparse storage to full dense storage.
 * intended only for verification.
 * ----------------------------------------------------------
 */
 const double dzero = 0;
 int ia = 0;
 int ja = 0;
 
 for(ja=0; ja < ncol_A; ja++) {
 for(ia=0; ia < nrow_A; ia++) {
   A(ia,ja) = dzero;
   };
   };

 for(ia=0; ia < nrow_A; ia++) {
   int istart = arowptr[ia];
   int iend = arowptr[ia+1]-1;
   int k = 0;
   for(k=istart; k <= iend; k++) {
     int ja = acol[k];
     double aij = aval[k];

     assert((0 <= ja) && (ja < ncol_A));

     A(ia,ja) = aij;
     };
   };

}
#undef A
