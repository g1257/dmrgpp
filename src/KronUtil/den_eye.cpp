#include "util.h"
void den_eye(       const int nrow_C, 
                    const int ncol_C, 

                    PsimagLite::Matrix<double>& c_ )
{
/*
 * -------------------------------------
 * form a dense identity matrix
 * -------------------------------------
 */

 int ic = 0;
 int jc = 0;

 for(jc=0; jc < ncol_C; jc++) {
 for(ic=0; ic < nrow_C; ic++) {

    c_(ic,jc) = (ic == jc) ? 1 : 0;
    };
    };

}

