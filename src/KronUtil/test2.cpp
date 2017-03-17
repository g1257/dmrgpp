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
     PsimagLite::CrsMatrix<double> a(a_);

	 PsimagLite::CrsMatrix<double> b(b_);

    /*
     * -----------------------------
     * explicitly form C = kron(A,B)
     * -----------------------------
     */

    const int nrow_C = nrow_A * nrow_B;
    const int ncol_C = ncol_A * ncol_B;
	if (ncol_C < 2) continue;
    PsimagLite::Matrix<double> c_(nrow_C, ncol_C);

    den_kron_form( nrow_A, ncol_A, a_,
                   nrow_B, ncol_B, b_,
                   c_);

    /*
     * ---------------------------------------
     * generate compressed sparse version of C
     * ---------------------------------------
     */
    PsimagLite::CrsMatrix<double> c(c_);


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

     const int max_nnz_D = 1 + den_nnz(d_);

	 PsimagLite::CrsMatrix<double> sd(d_.n_row(),d_.n_col());

     csr_submatrix( c,
                 
                    nrindex, ncindex, 
                    max_nnz_D,
                    rindex, cindex,
                    
                    sd );


     /*
      * -----------------------
      * convert to dense matrix
      * -----------------------
      */
      PsimagLite::Matrix<double> dd_(nrow_D, ncol_D);
	  crsMatrixToFullMatrix(dd_, sd);

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

    PsimagLite::CrsMatrix<double> e(nrow_E, ncol_E);
	const int max_nnz_E = 1 + den_nnz(e_);
    csr_kron_submatrix(
                a,
                b,
                nrindex,ncindex,
                max_nnz_E,
                rindex, cindex,
                e);

    /*
     * ---------------------------------
     * convert from sparse back to dense
     * ---------------------------------
     */
    PsimagLite::Matrix<double> se_(nrow_E, ncol_E);

	crsMatrixToFullMatrix(se_, e);


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
