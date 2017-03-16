#include "util.h"
#include "KronUtil.h"

int main()
{

  int nerrors = 0;
  double threshold = 0;
  int nrow_A = 0;
  int ncol_A = 0;
  int nrow_B = 0;
  int ncol_B = 0;

  for(threshold=0.01; threshold <= 1.1; threshold += 0.1) {
  for(ncol_A=1; ncol_A <= 10; ncol_A += 3 ) {
  for(nrow_A=1; nrow_A <= 10; nrow_A += 3 ) {
  for(ncol_B=1; ncol_B <= 10; ncol_B += 3 ) {
  for(nrow_B=1; nrow_B <= 10; nrow_B += 3 ) {

     PsimagLite::Matrix<double> a_(nrow_A, ncol_A);
     PsimagLite::Matrix<double> b_(nrow_B, ncol_B);

     den_gen_matrix( nrow_A, ncol_A,  threshold, a_ );
     den_gen_matrix( nrow_B, ncol_B,  threshold, b_ );

    /*
     * -----------------------------------
     * generate compressed row format from
     * dense matrices A, B
     * -----------------------------------
     */
     const int max_nnz_A =  1 + den_nnz(nrow_A,ncol_A, a_);
     PsimagLite::Vector<int>::Type arowptr(nrow_A+1);
     PsimagLite::Vector<int>::Type acol(max_nnz_A);
     PsimagLite::Vector<double>::Type aval(max_nnz_A);

     den2csr( nrow_A, ncol_A, a_,
              max_nnz_A, 
              arowptr, acol, aval );



     const int max_nnz_B = 1 + den_nnz(nrow_B,ncol_B, b_);
     PsimagLite::Vector<int>::Type browptr(nrow_B+1);
     PsimagLite::Vector<int>::Type bcol(max_nnz_B);
     PsimagLite::Vector<double>::Type bval(max_nnz_B);

     den2csr( nrow_B, ncol_B, b_,
              max_nnz_B, 
              browptr, bcol, bval );






    /*
     * -----------------------------
     * explicitly form C = kron(A,B)
     * -----------------------------
     */

    const int nrow_C = nrow_A * nrow_B;
    const int ncol_C = ncol_A * ncol_B;
    PsimagLite::Matrix<double> c_(nrow_C, ncol_C);

    den_kron_form( nrow_A, ncol_A, a_,
                   nrow_B, ncol_B, b_,
                   c_);

    /*
     * ---------------------------------------
     * generate compressed sparse version of C
     * ---------------------------------------
     */
    const int max_nnz_C = 1 + den_nnz( nrow_C, ncol_C, c_);
    PsimagLite::Vector<int>::Type scrowptr(nrow_C+1);
    PsimagLite::Vector<int>::Type sccol(max_nnz_C);
    PsimagLite::Vector<double>::Type scval(max_nnz_C);

    den2csr( nrow_C, ncol_C, c_,
             max_nnz_C,
             scrowptr, sccol,scval );



    PsimagLite::Vector<int>::Type rindex(nrow_C);
    PsimagLite::Vector<int>::Type cindex(ncol_C);
    
    int nrindex = 0;
    int ncindex = 0;

    /*
     * --------------------
     * extract  even rows
     * extract odd columns
     * --------------------
     */
    int ic = 0;
    int jc = 0;

    nrindex = 0;
    for(ic=0; ic < nrow_C; ic +=2 ) {
       rindex[nrindex] = ic;
       nrindex++;
       };

    ncindex = 0;
    for(jc=1; jc < ncol_C; jc += 2) {
       cindex[ncindex] = jc;
       ncindex++;
       };

    /*
     * -------------------------------
     * extract submatrix from C into D
     * -------------------------------
     */

    int nrow_D = nrindex;
    int ncol_D = ncindex;
    PsimagLite::Matrix<double> d_(nrow_D, ncol_D);

    den_submatrix( nrow_C, ncol_C, c_,
                   nrindex, ncindex,
                   rindex,  cindex,
                   d_);
    /*
     * -------------------------------------
     * generate submatrix from sparse version C
     * -------------------------------------
     */

     const int max_nnz_D = 1 + den_nnz( nrow_D,ncol_D, d_);
     PsimagLite::Vector<int>::Type sdrowptr(nrow_D + 1);
     PsimagLite::Vector<int>::Type sdcol(max_nnz_D);
     PsimagLite::Vector<double>::Type sdval(max_nnz_D);

     csr_submatrix( nrow_C, ncol_C, 
                    scrowptr,sccol,scval,
                 
                    nrindex, ncindex, 
                    max_nnz_D,
                    rindex, cindex,
                    
                    sdrowptr,sdcol,sdval );


     /*
      * -----------------------
      * convert to dense matrix
      * -----------------------
      */
      PsimagLite::Matrix<double> dd_(nrow_D, ncol_D);
      csr2den( nrow_D, ncol_D,
               sdrowptr, sdcol, sdval,
               dd_);

    /*
     * --------------------------------------
     * check both matrices should be the same
     * --------------------------------------
     */
     { 
     int id = 0;
     int jd = 0;

     for(jd=0; jd < ncol_D; jd++) {
     for(id=0; id < nrow_D; id++) {
         double diff = ABS( dd_(id,jd) - d_(id,jd) );
         const double tol = 1.0/(1000.0 * 1000.0);
         int isok = (diff <= tol );
         if (!isok) {
           nerrors += 1;
           printf("nrow_D %d ncol_D %d DD(%d,%d) %f D(%d,%d) %f diff %f\n",
                   nrow_D,ncol_D,   id,jd,dd_(id,jd),  id,jd,d_(id,jd), diff );

           };
        };
        };
     }
     
    
    /*
     * --------------------------------
     * form E = C(rindex(:),cindex(:))
     * without forming C
     * E should be the same as matrix D
     * --------------------------------
     */
    int nrow_E = nrindex;
    int ncol_E = ncindex;
    PsimagLite::Matrix<double> e_(nrow_E, ncol_E );

    den_kron_submatrix( nrow_A,ncol_A, a_,
                        nrow_B,ncol_B, b_,
                        nrindex, ncindex,
                        rindex,  cindex,
                        e_);

    /*
     * --------------------------
     * check E and D are the same
     * --------------------------
     */
    int ie = 0;
    int je = 0;
    for(je=0; je < ncol_E; je++) {
    for(ie=0; ie < nrow_E; ie++) {
       double eij = e_(ie,je);
       double dij = d_(ie,je);

       int isok = (eij == dij);
       if (!isok) {
         nerrors += 1;

         printf("nrow_A %d ncol_A %d nrow_B %d ncol_B %d nrindex %d ncindex %d\n",
                 nrow_A,  ncol_A, nrow_B, ncol_B, nrindex, ncindex );

         printf("ie %d je %d eij %f dij %f \n",
                 ie,je,  eij, dij );
         };
       };
       };



    /*
     * ------------------------------
     * form sparse version of E from 
     * sparse version of A, B
     * ------------------------------
     */
    const int max_nnz_E = 1 + den_nnz(nrow_E,ncol_E, e_);
    PsimagLite::Vector<int>::Type erowptr(nrow_E+1);
    PsimagLite::Vector<int>::Type ecol(max_nnz_E);
    PsimagLite::Vector<double>::Type eval(max_nnz_E);

    csr_kron_submatrix(
                nrow_A,ncol_A,arowptr,acol,aval,
                nrow_B,ncol_B,browptr,bcol,bval,
                nrindex,ncindex,
                max_nnz_E,
                rindex, cindex,
                erowptr,ecol,eval );

    /*
     * ---------------------------------
     * convert from sparse back to dense
     * ---------------------------------
     */
    PsimagLite::Matrix<double> se_(nrow_E, ncol_E);

    csr2den( nrow_E,ncol_E, erowptr, ecol,eval,
	     se_);


    /*
     * -----------------------------------------
     * check E(ie,je) and SE(ie,je) are the same
     * -----------------------------------------
     */
    {
    int ie = 0;
    int je = 0;

    for(je=0; je < ncol_E; je++) {
    for(ie=0; ie < nrow_E; ie++) {
      double diff = ABS( e_(ie,je) - se_(ie,je) );
      const double tol = 1.0/(1000.0 * 1000.0);
      int isok = (diff <= tol);
      if (!isok) {
        nerrors += 1;
        printf("nrow_E %d ncol_E %d E(%d,%d) %f SE(%d,%d) %f\n",
                nrow_E, ncol_E,  ie,je, e_(ie,je), ie,je, se_(ie,je) );

        };
     };
     };

    };

   };
   };
   };
   };
   };

 if (nerrors == 0) {
    printf("pass all tests\n");
    };
 return(0);
}

#undef X1
#undef X2
#undef X3
#undef SX1
#undef SX2
#undef SX3
#undef A
#undef B
#undef Y
