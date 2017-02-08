#include "util.h"

void csr_submatrix( const int nrow_A, 
                    const int ncol_A, 
                    const int arowptr[], 
                    const int acol[], 
                    const double aval[], 

                    const int nrow_B, 
                    const int ncol_B, 
                    const int max_nnz,

                    const int rindex[], 
                    const int cindex[],

                    int browptr[], 
                    int bcol[], 
                    double bval[] )
{
/*
 * ---------------------------------------------------------------------------
 * Input:   sparse matrix A in compressed sparse row format
 *
 * Input:   list of row index rindex(0:(nrow_B-1)) and
 *          list of column index cindex(0:(ncol_B-1))
 *
 * Out:     extract B = A( rindex(:), cindex(:)) in compress sparse row format
 * ---------------------------------------------------------------------------
 */

  
 /* 
  * ----------------------------------------------------
  * setup boolean array for fast mapping of column index
  * ----------------------------------------------------
  */
  int cmap[ncol_A];
  int nnz[nrow_B];

  int ja = 0;
  int ib = 0;
  int jb = 0;

  for(ja = 0; ja < ncol_A;  ja++) {
     cmap[ja] = -1;
     };

  for(jb=0; jb < ncol_B; jb++) {
     int ja = cindex[jb];
     assert((0 <= ja) && (ja < ncol_A));

     cmap[ ja ] = jb;
     };

  /*
   * -------------------------------------------------------
   * first pass to calculate number of nonzeros per row in B
   * -------------------------------------------------------
   */

  for(ib=0; ib < nrow_B; ib++) {
     nnz[ib] = 0;
     };


  for(ib=0; ib < nrow_B; ib++) {
     int ia = rindex[ib];
     int istart = arowptr[ia];
     int iend = arowptr[ia+1]-1;

     assert((0 <= ia) && (ia < nrow_A));

     int k = 0;
     for( k=istart; k <= iend; k++) {
         int ja = acol[k];
         assert((0 <= ja) && (ja < ncol_A));
         
         int jb = cmap[ ja ];
         int isok_jb = (0 <= jb) && (jb < ncol_B);
         if (isok_jb) {
             nnz[ib] += 1;
             };
         };
     };


  /*
   * ----------------------------
   * check for sufficient storage
   * ----------------------------
   */
  int total_nnz = 0;
  for(ib=0; ib < nrow_B; ib++) {
    total_nnz += nnz[ib];
    };
  assert( total_nnz <= max_nnz);

 /*
  * ---------------------------------
  * prefix sum to setup row pointers
  * ---------------------------------
  */
  
  browptr[0] = 0;
  for(ib=0; ib < nrow_B; ib++) {
     browptr[ib+1] = browptr[ib] + nnz[ib];
     };


  /*
   * ------------------------
   * reset array for 2nd pass
   * ------------------------
   */
  for(ib=0; ib < nrow_B; ib++) {
     nnz[ib] = 0;
     };


  /*
   * --------------------------------------
   * second pass to fill in compress sparse row 
   * data structure
   * --------------------------------------
   */
  for(ib=0; ib < nrow_B; ib++) {
      int ia = rindex[ib];
      int istart = arowptr[ia];
      int iend = arowptr[ia+1]-1;
      int k = 0;
      for(k=istart; k <= iend; k++) {
         int ja = acol[k];
         double aij = aval[k];

         jb = cmap[ja];
         int isvalid =  (0 <= jb) && (jb < ncol_B);
         if (isvalid) {
           int ipos = browptr[ib]  + nnz[ib];

           bval[ipos] = aij;
           bcol[ipos] = jb;
           nnz[ib] += 1;
           };
        };
      };

}
