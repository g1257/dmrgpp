#include "util.h"
void den_submatrix( const int nrow_A, 
                    const int ncol_A, 
                    const double a_[],

                    const int nrindex, 
                    const int ncindex,
                    const int rindex[],  
                    const int cindex[],

                    double c_[] )
{
/*
 * -------------------------------------
 * extract submatrix from a dense matrix
 * -------------------------------------
 */
 const int nrow_C = nrindex;
 const int ncol_C = ncindex;
#define C(ic,jc) c_[ (ic) + (jc)*nrow_C ]
#define A(ia,ja) a_[ (ia) + (ja)*nrow_A ]

 int ic = 0;
 int jc = 0;
 
 /*
  * ------------------------
  * check rindex[], cindex[]
  * ------------------------
  */
#ifndef NDEBUG
 for(ic=0; ic < nrindex; ic++) {
    int ia = rindex[ic];
    int isok = (0 <= ia) && (ia < nrow_A);
    assert( isok );
    };

 for(jc=0; jc < ncindex; jc++) {
    int ja = cindex[jc];
    int isok = (0 <= ja) && (ja < ncol_A);
    assert( isok );
    };
#endif

 for(jc=0; jc < ncol_C; jc++) {
 for(ic=0; ic < nrow_C; ic++) {
    int ia = rindex[ic];
    int ja = cindex[jc];

    C(ic,jc) = A(ia,ja);
    };
    };

}
 
#undef A
#undef C
