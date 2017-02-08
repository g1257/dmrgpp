#include "util.h"

void coord2csr( 
              const int nrow_A, 
              const int ncol_A, 
              const int nnz, 
              const int ilist[], 
              const int jlist[], 
              const double alist[],
              int arowptr[], 
              int acol[], 
              double aval[] )
{
/*
 * --------------------------------------
 * convert from coordinate storage format
 * "i j aij"
 * to compressed sparse row format
 *
 * need arowptr[] to have size nrow_A + 1
 * need acol[] to have size nnz
 * need aval[] to have size nnz
 * --------------------------------------
 */

 int nnz_row[nrow_A];


 {
 int ia = 0;
 for(ia=0; ia < nrow_A; ia++) {
    nnz_row[ia] = 0;
    };
 }

 /*
  * ------------------------
  * check ilist(:), jlist(:)
  * ------------------------
  */
#ifndef NDEBUG
 {
  int k = 0;
  for(k = 0; k < nnz; k++) {
      int ia = ilist[k];
      int ja = jlist[k];

     
      int isok_ia = (0 <= ia) && (ia < nrow_A);
      int isok_ja = (0 <= ja) && (ja < ncol_A);
       
      assert( isok_ia );
      assert( isok_ja );
      };
 }
#endif
      

 /*
  * -------------------------------------------------
  * first pass to compute number of non-zeros per row
  * -------------------------------------------------
  */
  {
  int k = 0;
  for(k=0; k < nnz; k++) {
    int ia = ilist[k];
    nnz_row[ia] += 1;
    };
  }

 /*
  * ------------------------------
  * prefix sum to setup row poiner 
  * ------------------------------
  */
  {
  int ia = 0;

  arowptr[0] = 0;
  for(ia=0; ia < nrow_A; ia++) {
      arowptr[ia+1] = arowptr[ia] + nnz_row[ia];
      };

  for(ia=0; ia < nrow_A; ia++) {
    nnz_row[ia] = 0;
    };

  }

 /*
  * -------------------------------------------------------------
  * second pass to fill in data into compressed sparse row format
  * -------------------------------------------------------------
  */
  {
  int k = 0;

  for(k=0; k < nnz; k++) {
    int ia = ilist[k];
    int ja = jlist[k];
    double aij = alist[k];

    int ipos = arowptr[ia] + nnz_row[ia];

    acol[ipos] = ja;
    aval[ipos] = aij;

    nnz_row[ia] += 1;
    };
  }

}
