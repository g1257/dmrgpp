#include "util.h"

void den_copymat( const int nrow, 
                  const int ncol, 
                  const int asrc_[], 
                        int bdest_[] )
{
/*
 * -----------------------------
 * copy  matrix B(:,:) =  A(:,:)
 * -----------------------------
 */

#define Asrc(ia,ja) asrc_[ (ia) + (ja)*nrow ]
#define Bdest(ib,jb) bdest_[ (ib) + (jb)*nrow ]

  int ia = 0;
  int ja = 0;
  
  for(ja=0; ja < ncol; ja++) {
  for(ia=0; ia < nrow; ia++) {
     Bdest(ia,ja) = Asrc(ia,ja);
     };
     };
}
#undef Asrc
#undef Bdest

