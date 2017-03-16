#include "util.h"
void csr2den(  const int nrow_A,
               const int ncol_A, 
               const PsimagLite::Vector<int>::Type& arowptr,
               const PsimagLite::Vector<int>::Type& acol,
               const PsimagLite::Vector<double>::Type& aval,

               PsimagLite::Matrix<double>& a_)
{
/*
 * ----------------------------------------------------------
 * convert from compress sparse storage to full dense storage.
 * intended only for verification.
 * ----------------------------------------------------------
 */
 a_.setTo(0.0);

 for(int ia=0; ia < nrow_A; ia++) {
   int istart = arowptr[ia];
   int iend = arowptr[ia+1]-1;
   int k = 0;
   for(k=istart; k <= iend; k++) {
     int ja = acol[k];
     double aij = aval[k];

     assert((0 <= ja) && (ja < ncol_A));

     a_(ia,ja) = aij;
     };
   };

}

