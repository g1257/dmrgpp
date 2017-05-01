#include "util.h"
bool den_is_eye(const PsimagLite::Matrix<double>& a_)
{
	const int nrow_A = a_.n_row();
	const int ncol_A = a_.n_col();
/*
 * -------------------------
 * return whether A is the identity matrix
 * matrix A in dense storage format
 * -------------------------
 */


  if (nrow_A != ncol_A) { return( false ); };


  int ja = 0;
  for(ja=0; ja < ncol_A; ja++) {
    int ia = 0;
    for(ia=0; ia < nrow_A; ia++) {
       double aij = a_(ia,ja);
       double eij = (ia == ja) ? 1 : 0;

       bool is_eye = (aij == eij);
       if (!is_eye) { return( false ); };
       };
    };

   /*
    * -----------------
    * passed all checks
    * -----------------
    */
   return( true );

}
#undef A
