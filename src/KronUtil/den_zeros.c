#include "util.h"

void den_zeros( const int nrow_A, 
                const int ncol_A, 
                      double a_[] )
{
#define A(ia,ja) a_[(ia) + (ja)*nrow_A]
/*
 * ----------------------
 * set array to all zeros
 * ----------------------
 */
 int ia = 0;
 int ja = 0;
 const double dzero = 0;

 for(ja=0; ja < ncol_A; ja++) {
 for(ia=0; ia < nrow_A; ia++) {
   A(ia,ja) = dzero;
   };
   };
}
#undef A
