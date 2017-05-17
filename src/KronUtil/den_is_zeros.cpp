#include "util.h"

template<typename ComplexOrRealType>
bool den_is_zeros(const PsimagLite::Matrix<ComplexOrRealType>& a_)
{
	const int nrow_A = a_.n_row();
	const int ncol_A = a_.n_col();
/*
 * -------------------------
 * return whether A is the zero matrix
 * matrix A in dense storage format
 * -------------------------
 */

  const ComplexOrRealType zero = 0;



  int ja = 0;
  for(ja=0; ja < ncol_A; ja++) {
    int ia = 0;
    for(ia=0; ia < nrow_A; ia++) {
       ComplexOrRealType aij = a_(ia,ja);
       if (aij != zero) { return( false ); };

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
