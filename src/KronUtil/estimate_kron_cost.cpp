#include "util.h"



void estimate_kron_cost( const int nrow_A,
                         const int ncol_A,
                         const int nnz_A_in, 
                         const int nrow_B,
                         const int ncol_B,
                         const int nnz_B_in,
                         double *p_kron_nnz, 
                         double *p_kron_flops, 
                         int *p_imethod )
{
  int imethod = 1;
  
  assert( p_kron_nnz != 0 );
  assert( p_kron_flops != 0 );
  assert( p_imethod != 0 );
  
  double nnz_A = (double) nnz_A_in;
  double nnz_B = (double) nnz_B_in;
  bool is_dense_A = (nnz_A_in == nrow_A * ncol_A);
  bool is_dense_B = (nnz_B_in == nrow_B * ncol_B);


  /*
   * ------------------------------------
   * assume dense matrix operations are
   * faster than sparse matrix operations
   * ------------------------------------
   */
  const double dense_flop_discount = 0.2;
  double discount = 1;
  
  /*
  % ------------------------
  % method 1:
  %
  % BY(ib,ja) = B(ib,jb)*Y( jb,ja ),
  %
  % BY has size  nrow_B by ncol_A
  %
  % X(ib,ia) = BY(ib,ja)*At(ja,ia)
  % ------------------------
  */
  discount = (is_dense_B) ? dense_flop_discount : 1;
  double flops_BY =  discount * 2.0*nnz_B*ncol_A;

  discount = (is_dense_A) ? dense_flop_discount : 1;
  double flops_BYAt = discount * 2.0*nnz_A*nrow_B;

  double flops_method1 = flops_BY + flops_BYAt;
  double nnz_method1 = nrow_B * ncol_A;


  /*
  % ------------------------
  % method 2:
  %
  % YAt(jb,ia)  = Y(jb,ja) * At(ja,ia)
  %
  % YAt has size ncol_B by nrow_A
  %
  % X(ib,ia) = B(ib,jb)*YAt(jb,ia)
  % ------------------------
  */
  discount = (is_dense_A) ? dense_flop_discount : 1;
  double flops_YAt = discount * 2.0*nnz_A*ncol_B;

  discount = (is_dense_B) ? dense_flop_discount : 1;
  flops_BYAt = discount * 2.0*nnz_B*nrow_A;

  double flops_method2 = flops_YAt + flops_BYAt;
  double nnz_method2 = ncol_B * nrow_A;

  
  /*
  % ------------------
  % method 3:
  % X = kron(A,B) * Y
  %
  % evaluate all nonzeros in A and nonzeros in B
  % need A and B to be very sparse
  % ------------------
  */
  double flops_method3 =  2.0*nnz_A*nnz_B;
  double nnz_method3 = nnz_A*nnz_B;

  double kron_flops = MIN( flops_method1, 
                            MIN( flops_method2, flops_method3));

  double kron_nnz = MIN( nnz_method1, 
                          MIN( nnz_method2, nnz_method3));

  const int minimize_flops = TRUE;
  if (minimize_flops) {
    /*
     * ----------------------------
     * minimize the number of flops
     * ----------------------------
     */

    imethod = 3;
    if (kron_flops == flops_method1) {
        imethod = 1;
        };
    if (kron_flops == flops_method2) {
        imethod = 2;
        };
   }
  else {
    /*
     * ------------------------------------
     * minimize amount of temporary storage
     * ------------------------------------
     */
    imethod = 3;
    if (kron_nnz == nnz_method1) {
          imethod = 1;
          };
    if (kron_nnz == nnz_method2) {
          imethod = 2;
          };
    }

  

  
  
  *p_kron_nnz = kron_nnz;
  *p_kron_flops = kron_flops;
  *p_imethod = imethod;
}


#undef MIN
#undef TRUE
#undef FALSE
