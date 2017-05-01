#include "util.h"
int den_nnz(const PsimagLite::Matrix<double>& a_)
{
	const int nrow_A = a_.n_row();
	const int ncol_A = a_.n_col();
/*
 * -------------------------
 * return number of nonzeros
 * matrix A in dense storage format
 * -------------------------
 */
  const bool use_estimate = true;
  if (use_estimate) {
     return( nrow_A * ncol_A );
     };

  

  const double dzero = 0;
  int nnz_A = 0;

  int ja = 0;

  /*
   * -----------------------------
   * Note that  nnz_A is a reduction variable
   * -----------------------------
   */
  for(ja=0; ja < ncol_A; ja++) {
    int ia = 0;
    for(ia=0; ia < nrow_A; ia++) {
       int is_zero = (a_(ia,ja) == dzero);
       nnz_A = (is_zero)? nnz_A : (nnz_A+1);
       };
     };

  return( nnz_A );
}
#undef A
