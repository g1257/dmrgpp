#include "util.h"
#include "KronUtil.h"

#ifndef USE_FLOAT
typedef double RealType;
#else
typedef float RealType;
#endif

int main()
{
  typedef RealType RealType;
  typedef std::complex<RealType> ComplexOrRealType;
  int nerrors = 0;
  RealType thresholdA = 0;
  RealType thresholdB = 0;
  int nrow_A = 0;
  int ncol_A = 0;
  int nrow_B = 0;
  int ncol_B = 0;

  for(thresholdB=0; thresholdB <= 1.1; thresholdB += 0.1) {
  for(thresholdA=0; thresholdA <= 1.1; thresholdA += 0.1) {
  for(ncol_A=1; ncol_A <= 10; ncol_A += 3 ) {
  for(nrow_A=1; nrow_A <= 10; nrow_A += 3 ) {
  for(ncol_B=1; ncol_B <= 10; ncol_B += 3 ) {
  for(nrow_B=1; nrow_B <= 10; nrow_B += 3 ) {

     PsimagLite::Matrix<ComplexOrRealType> a_(nrow_A, ncol_A);
     PsimagLite::Matrix<ComplexOrRealType> b_(nrow_B, ncol_B);


     if ((thresholdA == 0) && (nrow_A == ncol_A)) {
       /*
        * ------------------------------------
        * special case to test identity matrix
        * ------------------------------------
        */
       den_eye( nrow_A, ncol_A, a_ );
       assert( den_is_eye(a_) );

       PsimagLite::CrsMatrix<ComplexOrRealType> a(a_);
       assert( csr_is_eye(a) );
       }
     else {
       den_gen_matrix( nrow_A, ncol_A,  thresholdA, a_);

       PsimagLite::CrsMatrix<ComplexOrRealType> a(a_);
       assert( den_is_eye(a_) == csr_is_eye(a) );
       assert( den_is_zeros(a_) == csr_is_zeros(a) );

       };


     if ((thresholdB == 0) && (nrow_B == ncol_B)) {
       den_eye( nrow_B, ncol_B, b_ );
       assert( den_is_eye(b_));

       PsimagLite::CrsMatrix<ComplexOrRealType> b(b_);
       assert( csr_is_eye(b) );
       }
     else {
       den_gen_matrix( nrow_B, ncol_B,  thresholdB, b_);

       PsimagLite::CrsMatrix<ComplexOrRealType> b(b_);
       assert( den_is_eye(b_) == csr_is_eye(b) );
       assert( den_is_zeros(b_) == csr_is_zeros(b) );
       };
    /*
     * -----------------------------------
     * generate compressed row format from
     * dense matrices A, B
     * -----------------------------------
     */
     PsimagLite::CrsMatrix<ComplexOrRealType> a(a_);
     assert( den_is_eye(a_) == csr_is_eye(a) );

     PsimagLite::CrsMatrix<ComplexOrRealType> b(b_);
     assert( den_is_eye(b_) == csr_is_eye(b) );

    /*
     * -----------------------------
     * explicitly form C = kron(A,B)
     * -----------------------------
     */

    const int nrow_C = nrow_A * nrow_B;
    const int ncol_C = ncol_A * ncol_B;
	if (ncol_C < 2) continue;
    PsimagLite::Matrix<ComplexOrRealType> c_(nrow_C, ncol_C);

    den_kron_form( nrow_A, ncol_A, a_,
                   nrow_B, ncol_B, b_,
                   c_);

    /*
     * ---------------------------------------
     * generate compressed sparse version of C
     * ---------------------------------------
     */
    PsimagLite::CrsMatrix<ComplexOrRealType> c(c_);


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
    PsimagLite::Matrix<ComplexOrRealType> d_(nrow_D, ncol_D);

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

	 PsimagLite::CrsMatrix<ComplexOrRealType> sd(d_.n_row(),d_.n_col());

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
      PsimagLite::Matrix<ComplexOrRealType> dd_(nrow_D, ncol_D);
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
         RealType diff = std::abs(dd_(id,jd) - d_(id,jd) );
         const RealType tol = 1.0/(1000.0 * 1000.0);
         bool isok = (diff <= tol );
         if (!isok) {
           nerrors += 1;
           std::cout<<"nrow_D "<<nrow_D;
		   std::cout<<" ncol_D "<<ncol_D;
		   std::cout<<" DD("<<id<<","<<jd<<")"<<dd_(id, jd);
		   std::cout<<" D("<<id<<","<<jd<<")"<<d_(id, jd);
		   std::cout<<" diff "<<diff<<"\n";
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
    PsimagLite::Matrix<ComplexOrRealType> e_(nrow_E, ncol_E );

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
       ComplexOrRealType eij = e_(ie,je);
       ComplexOrRealType dij = d_(ie,je);

       int isok = (eij == dij);
       if (!isok) {
         nerrors += 1;

         std::cout<<"nrow_A "<<nrow_A<<" ncol_A "<<ncol_A;
		 std::cout<<" nrow_B "<<nrow_B<<" ncol_B "<<ncol_B;
		 std::cout<<"  nrindex "<<nrindex<<" ncindex "<<ncindex<<"\n";

		 std::cout<<"ie "<<ie<< " je "<<je<<" eij "<<eij<<" dij "<<dij<<"\n";
         }
       }
       }


    /*
     * ------------------------------
     * form sparse version of E from
     * sparse version of A, B
     * ------------------------------
     */

    PsimagLite::CrsMatrix<ComplexOrRealType> e(nrow_E, ncol_E);
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
    PsimagLite::Matrix<ComplexOrRealType> se_(nrow_E, ncol_E);

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
      RealType diff = std::abs( e_(ie,je) - se_(ie,je) );
      const RealType tol = 1.0/(1000.0 * 1000.0);
      bool isok = (diff <= tol);
      if (!isok) {
        nerrors += 1;
        std::cout<<"nrow_E "<<nrow_E<<" ncol_E "<<ncol_E;
		std::cout<<" E("<<ie<<","<<je<<" "<<e_(ie, je);
		std::cout<<" SE("<<ie<<","<<je<<")"<<se_(ie, je)<<"\n";
        }
     }
     }

    }

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
