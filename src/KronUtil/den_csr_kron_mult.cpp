#include "util.h"

template<typename ComplexOrRealType>
void den_csr_kron_mult_method(const int imethod,
                              const char transA,
                              const char transB,
                              const PsimagLite::Matrix<ComplexOrRealType>& a_,
                              const PsimagLite::CrsMatrix<ComplexOrRealType>& b,
                              const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                              SizeType offsetY,
                              typename PsimagLite::Vector<ComplexOrRealType>::Type& xout_,
                              SizeType offsetX)
{
	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;

	const int isTransA = (transA == 'T') || (transA == 't');
	const int isTransB = (transB == 'T') || (transB == 't');

	const int isConjTransA = (transA == 'C') || (transA == 'c');
	const int isConjTransB = (transB == 'C') || (transB == 'c');

	const int nrow_A = a_.rows();
	const int ncol_A = a_.cols();
	const int nrow_B = b.rows();
	const int ncol_B = b.cols();

	const int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	const int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;
	const int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	const int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

	const int nrow_X = nrow_2;
	const int ncol_X = nrow_1;
	const int nrow_Y = ncol_2;
	const int ncol_Y = ncol_1;

	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_Y, ncol_Y, yin_, offsetY);
	PsimagLite::MatrixNonOwned<ComplexOrRealType> xout(nrow_X, ncol_X, xout_, offsetX);

	assert((imethod == 1) ||
	       (imethod == 2) ||
	       (imethod == 3));

	bool no_work = (den_is_zeros(a_) || csr_is_zeros(b));
	if (no_work) {
		return;
	};
/*
 *   -------------------------------------------------------------
 *   A in dense matrix storage format
 *   B in compressed sparse ROW format
 *
 *   X += kron( op(A), op(B)) * Y
 *   X +=  op(B) * Y * transpose(op(A))
 *
 *   nrow_X = nrow_2,   ncol_X = nrow_1
 *   nrow_Y = ncol_2,   nrow_Y = nrow_2
 *
 *   that can be computed as either
 *   imethod == 1
 *
 *   X(ib,ia) +=  (B(ib,jb) * Y( jb,ja)) * transpose( A(ia,ja))
 *
 *   X(ix,jx) +=  (B2(ib2,jb2) * Y(iy,jy)) * transpose(A1(ia1,ja1))
 *
 *
 *
 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops
 *
 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
 *
 *   imethod == 2
 *
 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))
 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
 *
 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
 *
 *   imethod == 3
 *
 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B
 *
 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
 *   -------------------------------------------------------------
 */



	if (imethod == 1) {

	/*
	 *  --------------------------------------------
	 *  BY(ib,ja) = (B(ib,jb))*Y(jb,ja)
	 *
	 *  X(ib,ia) += BY(ib,ja ) * transpose(A(ia,ja))
	 *  --------------------------------------------
	 */

		int nrow_BY = nrow_X;
		int ncol_BY = ncol_Y;
		PsimagLite::Matrix<ComplexOrRealType> by_(nrow_BY, ncol_BY);
		PsimagLite::MatrixNonOwned<ComplexOrRealType> byRef(by_);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> byConstRef(by_);

	/*
	 * ---------------
	 * setup BY(nrow_BY,ncol_BY)
	 * ---------------
	 */

		{
			int iby = 0;
			int jby = 0;

			// not needed FIXME
			for(jby=0; jby < ncol_BY; jby++) {
				for(iby=0; iby < nrow_BY; iby++) {
					by_(iby,jby) = 0;
				};
			};
		}



		{
			/*
	 * ------------------------------
	 * BY(ib,ja)  = B(ib,jb)*Y(jb,ja)
	 * ------------------------------
	 */
			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;
			csr_matmul_pre(  trans,
			                 b,

			                 nrow_Y,
			                 ncol_Y,
			                 yin,

			                 nrow_BY,
			                 ncol_BY,
			                 byRef);

		}

		{
	/*
	 * -------------------------------------------
	 * X(ib,ia) += BY(ib,ja) * transpose(A(ia,ja))
	 * -------------------------------------------
	 */
			const char trans = (isTransA || isConjTransA) ? 'N' : 'T';
			if (is_complex && isConjTransA) {
				// -----------------------------------------
				// transpose( conj(transpose(A) ) is conj(A)
				// perform  conj operation
				// -----------------------------------------
                              PsimagLite::Matrix<ComplexOrRealType> a_conj(nrow_A,ncol_A);

			      for(int ja=0; ja < ncol_A; ja++) {
			      for(int ia=0; ia < nrow_A; ia++) {
				     a_conj(ia,ja) = PsimagLite::conj( a_(ia,ja) );
			      };
			      };


			      den_matmul_post(
			            trans,
			            nrow_A,
			            ncol_A,
			            a_conj,

			            nrow_BY,
			            ncol_BY,
			            byConstRef,

			            nrow_X,
			            ncol_X,
			            xout);


			}
			else {

			     den_matmul_post(
			            trans,
			            nrow_A,
			            ncol_A,
			            a_,

			            nrow_BY,
			            ncol_BY,
			            byConstRef,

			            nrow_X,
			            ncol_X,
			            xout);
			};


		}
	}
	else if (imethod == 2) {
		/*
	 * ---------------------
	 * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja))
	 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
	 * ---------------------
	 */

		int nrow_YAt = nrow_Y;
		int ncol_YAt = ncol_X;
		PsimagLite::Matrix<ComplexOrRealType> yat_(nrow_YAt, ncol_YAt);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> yatConstRef(yat_);
		PsimagLite::MatrixNonOwned<ComplexOrRealType> yatRef(yat_);

		/*
	 * ----------------
	 * setup YAt(jb,ia)
	 * ----------------
	 */

		{
			int iy = 0;
			int jy = 0;

			for(jy=0; jy < ncol_YAt; jy++) {
				for(iy=0; iy < nrow_YAt; iy++) {
					yat_(iy,jy) = 0;
				};
			};
		}




		{
			/*
	  * ---------------------
	  * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja)
	  * ---------------------
	  */
			const char transa = (isTransA || isConjTransA) ? 'N' : 'T';
			if (is_complex && isConjTransA) {
				// -----------------------------------------
				// transpose( conj(transpose(A) ) is conj(A)
				// perform in-place conj operation
				// -----------------------------------------
                              PsimagLite::Matrix<ComplexOrRealType> a_conj(nrow_A,ncol_A);
			      for(int ja=0; ja < ncol_A; ja++) {
                              for(int ia=0; ia < nrow_A; ia++) {
				      a_conj(ia,ja) = PsimagLite::conj( a_(ia,ja) );
			      };
			      };

			      den_matmul_post( transa,
			                 nrow_A,
			                 ncol_A,
			                 a_conj,

			                 nrow_Y,
			                 ncol_Y,
			                 yin,

			                 nrow_YAt,
			                 ncol_YAt,
			                 yatRef);

			}
			else {

			      den_matmul_post( transa,
			                 nrow_A,
			                 ncol_A,
			                 a_,

			                 nrow_Y,
			                 ncol_Y,
			                 yin,

			                 nrow_YAt,
			                 ncol_YAt,
			                 yatRef);
			};


		}




		{
	/*
	 * ------------
	 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
	 * ------------
	 */

			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;
			csr_matmul_pre( trans,
			                b,

			                nrow_YAt,
			                ncol_YAt,
			                yatConstRef,

			                nrow_X,
			                ncol_X,
			                xout);

		}
	}
	else if (imethod == 3) {
	/*
	* ---------------------------------------------
	* C = kron(A,B)
	* C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
	* X([ib,ia]) += C([ib,ia],[jb,ja]) * Y([jb,ja])
	* ---------------------------------------------
	*/
		const bool use_inner_sparse = false;

		if (use_inner_sparse) {
	/*
	* -----------------------------
	* sparse matrix B in inner loop
	* -----------------------------
	*/
			int ia = 0;
			int ja = 0;
			for(ia=0; ia < nrow_A; ia++) {
				for(ja=0; ja < ncol_A; ja++) {
					ComplexOrRealType aij = a_(ia,ja);
					if (is_complex && isConjTransA) {
						aij = PsimagLite::conj( aij );
					};

					int ib = 0;
					for(ib=0; ib < nrow_B; ib++) {
						int istartb = b.getRowPtr(ib);
						int iendb = b.getRowPtr(ib + 1);

						int kb = 0;
						for(kb=istartb; kb < iendb; kb++) {
							int jb = b.getCol(kb);
							ComplexOrRealType bij = b.getValue(kb);
							if (is_complex && isConjTransB) {
								bij = PsimagLite::conj( bij );
							};

							ComplexOrRealType cij = aij * bij;

							int ix = (isTransB || isConjTransB) ? jb : ib;
							int jx = (isTransA || isConjTransA) ? ja : ia;
							int iy = (isTransB || isConjTransB) ? ib : jb;
							int jy = (isTransA || isConjTransA) ? ia : ja;

							xout(ix,jx) +=  (cij * yin(iy,jy));
						};
					};
				};
			};

		}
		else {
	/*
	* -----------------------------
	* sparse matrix B in outer loop
	* -----------------------------
	*/

			int ib = 0;
			for(ib=0; ib < nrow_B; ib++) {
				int istartb = b.getRowPtr(ib);
				int iendb = b.getRowPtr(ib + 1);

				int kb = 0;
				for(kb=istartb; kb < iendb; kb++) {
					int jb = b.getCol(kb);
					ComplexOrRealType bij = b.getValue(kb);
					if (is_complex && isConjTransB) {
						bij = PsimagLite::conj( bij );
					};

					int ia = 0;
					int ja = 0;

					for(ia=0; ia < nrow_A; ia++) {
						for(ja=0; ja < ncol_A; ja++) {

							ComplexOrRealType aij = a_(ia,ja);
							if (is_complex && isConjTransA) {
								aij = PsimagLite::conj( aij );
							};

							ComplexOrRealType cij = aij * bij;

							int ix = (isTransB || isConjTransB) ? jb : ib;
							int jx = (isTransA || isConjTransA) ? ja : ia;
							int iy = (isTransB || isConjTransB) ? ib : jb;
							int jy = (isTransA || isConjTransA) ? ia : ja;

							xout(ix,jx) +=  (cij * yin(iy,jy));
						};
					};
				};
			};

		};
	};

}

template<typename ComplexOrRealType>
void den_csr_kron_mult(const char transA,
                       const char transB,
                       const PsimagLite::Matrix<ComplexOrRealType>& a_,
                       const PsimagLite::CrsMatrix<ComplexOrRealType>& b,
                       const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                       SizeType offsetY,
                       typename PsimagLite::Vector<ComplexOrRealType>::Type& xout_,
                       SizeType offsetX,
                       const typename PsimagLite::Real<ComplexOrRealType>::Type denseFlopDiscount)

{
	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int idebug = 0;
/*
 *   -------------------------------------------------------------
 *   A in dense matrix format
 *   B in compressed sparse ROW format
 *
 *   X += kron( A, B) * Y
 *   that can be computed as either
 *   imethod == 1
 *
 *   X(ib,ia) += (B(ib,jb) * Y(jb,ja) ) * transpose(A(ia,ja)   or
 *               BY(ib,ja) = B(ib,jb)*Y(jb,ja)
 *               BY is nrow_B by ncol_A, need   2*nnz(B)*ncolA flops
 *
 *   X(ib,ia) +=   BY(ib,ja) * transpose(A(ia,ja)) need 2*nnz(A)*nrowB flops
 *
 *   imethod == 2
 *
 *   X(ib,ia) += B(ib,jb) * (Y(jb,ja) * transpose(A))    or
 *                YAt(jb,ia) = Y(jb,ja) * transpose(A(ia,ja))
 *                YAt is ncolB by nrowA, need 2*nnz(A) * ncolB flops
 *
 *   X(ib,ia) += B(ib,jb) * YAt(jb,ia)  need nnz(B) * nrowA flops
 *
 *   imethod == 3
 *
 *   X += kron(A,B) * Y   by visiting all non-zero entries in A, B
 *
 *   this is feasible only if A and B are very sparse, need nnz(A)*nnz(B) flops
 *   -------------------------------------------------------------
 */

	const int nrow_A = a_.rows();
	const int ncol_A = a_.cols();
	const int nrow_B = b.rows();
	const int ncol_B = b.cols();

	int nnz_A = den_nnz(a_);
	int nnz_B = csr_nnz(b);
	int no_work = (den_is_zeros(a_) || csr_is_zeros(b));
	if (no_work) {
		if (idebug >= 1) {
			printf("den_csr:no_work,nrow_A=%d,nrow_B=%d,nrow_B=%d,ncol_B=%d\n",
			       nrow_A,ncol_A,   nrow_B,ncol_B );
		};

		return;
	};

	const int isTransA = (transA == 'T') || (transA == 't');
	const int isTransB = (transB == 'T') || (transB == 't');

	const int isConjTransA = (transA == 'C') || (transA == 'c');
	const int isConjTransB = (transB == 'C') || (transB == 'c');

	/*
  * --------------------------
  * check for special case
  * that B is identity matrix
  * --------------------------
  */
	if (csr_is_eye(b)) {
		if (idebug >= 1) {
			printf("den_csr: B is eye: nrow_A=%d,ncol_A=%d,nrow_B=%d,ncol_B=%d\n",
			       nrow_A,ncol_A, nrow_B, ncol_B );

		};

	/*
	----------------------
	X +=  ( op(B) ) Y * transpose( op(A) )
	----------------------
	*/
		const int nrow_Y = (isTransB || isConjTransB) ? nrow_B : ncol_B;
		const int ncol_Y = (isTransA || isConjTransA) ? nrow_A : ncol_A;

		const int nrow_X = (isTransB || isConjTransB) ? ncol_B : nrow_B;
		const int ncol_X = (isTransA || isConjTransA) ? ncol_A : nrow_A;

		PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_Y, ncol_Y, yin_, offsetY);
		PsimagLite::MatrixNonOwned<ComplexOrRealType> xout(nrow_X, ncol_X, xout_, offsetX);



		const char  trans1 =  (isTransA || isConjTransA) ? 'N' : 'T';
		if (is_complex && isConjTransA) {
			// ---------------------------------------
			// transpose( conj(transpose(A)) is conj(A)
			// perform in-place conj operation
			// ---------------------------------------
                        PsimagLite::Matrix<ComplexOrRealType> a_conj(nrow_A,ncol_A);

			for(int ja=0; ja < ncol_A; ja++) {
			for(int ia=0; ia < nrow_A; ia++) {
				a_conj(ia,ja) = PsimagLite::conj( a_(ia,ja) );
			};
			};

		        den_matmul_post(  trans1,
		                  nrow_A, ncol_A, a_conj,
		                  nrow_Y, ncol_Y, yin,
		                  nrow_X, ncol_X, xout );
		}
		else {

		        den_matmul_post(  trans1,
		                  nrow_A, ncol_A, a_,
		                  nrow_Y, ncol_Y, yin,
		                  nrow_X, ncol_X, xout );
		};


		return;
	};

	ComplexOrRealType kron_nnz = 0;
	ComplexOrRealType kron_flops = 0;
	int imethod = 1;


	const int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	const int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;

	const int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	const int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;


	estimate_kron_cost( nrow_1,ncol_1,nnz_A, nrow_2,ncol_2,nnz_B,
	                    &kron_nnz, &kron_flops, &imethod, denseFlopDiscount);


	den_csr_kron_mult_method(
	            imethod,
	            transA,
	            transB,

	            a_,

	            b,

	            yin_,
	            offsetY,
	            xout_,
	            offsetX);
}

#undef A
#undef BY
#undef YAt
#undef X
#undef Y
