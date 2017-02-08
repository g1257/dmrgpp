#include "util.h"
void den_kron_form( const int nrow_A, 
                    const int ncol_A, 
                    const double a_[], 

                    const int nrow_B, 
                    const int ncol_B, 
                    const double b_[],   

                          double c_[] )
{
/*
 * ---------------------------------------
 * form C = kron(A,B),  where
 * nrow_C = nrow_A * nrow_B
 * ncol_C = ncol_A * ncol_B
 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
 * ---------------------------------------
 */
#define A(ia,ja) a_[(ia) + (ja)*nrow_A]
#define B(ib,jb) b_[(ib) + (jb)*nrow_B]
#define C(ib,jb) c_[(ic) + (jc)*nrow_C]

  const int nrow_C = nrow_A * nrow_B;

  int ia = 0;
  int ja = 0;
  int ib = 0;
  int jb = 0;


  for(ja=0; ja < ncol_A; ja++) {
  for(jb=0; jb < ncol_B; jb++) {
  for(ia=0; ia < nrow_A; ia++) {
  for(ib=0; ib < nrow_B; ib++) {
           int ic = ib + ia*nrow_B;
           int jc = jb + ja*ncol_B;

           C(ic,jc) = A(ia,ja) * B(ib,jb);
    };
    };
    };
    };
}
