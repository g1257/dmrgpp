#include "util.h"

void csr_transpose( const int nrow_A, 
                   const int ncol_A, 
                   const int arowptr[], 
                   const int acol[], 
                   const double aval[],  

                         int atrowptr[], 
                         int atcol[], 
                         double atval[] )
{
/*
 * --------------------------------------------------------
 * At = tranpose(A) where A in sparse compressed row format
 * --------------------------------------------------------
 */
 
 const int nrow_At = ncol_A;

 int nnz_row_At[nrow_At];
 
 {
 int iat = 0;
 for(iat=0; iat < nrow_At; iat++) {
   nnz_row_At[iat] = 0;
   };
 }

  
 /*
  * --------------------------------------
  * first pass to count number of nonzeros
  * per row in At = transpose(A)
  * --------------------------------------
  */

 { 
 int ia = 0;
 for(ia=0; ia < nrow_A; ia++ ) {
   int istart = arowptr[ia];
   int iend = arowptr[ia+1]-1;
   int k = 0;
   for(k=istart; k <= iend; k++) {
     int ja = acol[k];
     int iat = ja;

     nnz_row_At[iat] += 1;
     };
   };
  }

 /*
  * ---------------------------------------
  * prefix sum to setup row pointers for At
  * ---------------------------------------
  */
 {
 int iat = 0;
 atrowptr[0] = 0;
 for(iat=0; iat < nrow_At; iat++) {
   atrowptr[iat+1] = atrowptr[iat] + nnz_row_At[iat];
   };

 for(iat=0; iat < nrow_At; iat++) {
   nnz_row_At[iat] = 0;
   };
 }


 /*
  * ----------------------
  * second pass to fill At
  * ----------------------
  */
 {
 int ia = 0;
 for(ia=0; ia < nrow_A; ia++) {
   int istart = arowptr[ia];
   int iend = arowptr[ia+1]-1;
   int k = 0;
   for(k=istart; k <= iend; k++) {
    int ja = acol[k];
    double aij = aval[k];

    int iat = ja;
    int jat = ia;

    int ipos = atrowptr[iat] + nnz_row_At[iat];
    atcol[ipos] = jat;
    atval[ipos] = aij;

    nnz_row_At[iat] += 1;
    };
  };
 }

}


