#include "util.h"
void den_transpose( const int nrow_A, 
                    const int ncol_A, 
                    const double a_[],  
                          double at_[] )
{
/*
 * ---------------------------
 * perform transpose operation
 * for dense matrix A
 * ---------------------------
 */
#define A(ia,ja) a_[ (ia) + (ja)*nrow_A ]
#define At(ia,ja) at_[ (ia) + (ja)*nrow_At ]

const int nrow_At = ncol_A;

int ia = 0;
int ja = 0;

for(ja=0; ja < ncol_A; ja++) {
for(ia=0; ia < nrow_A; ia++) {
   At(ja,ia) = A(ia,ja);
   };
   };

}
#undef A
#undef At

