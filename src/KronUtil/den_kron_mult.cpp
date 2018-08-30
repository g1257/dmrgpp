#include "util.h"

template<typename ComplexOrRealType>
void den_kron_mult_method(const int imethod,
                          const char transA,
                          const char transB,
                          const PsimagLite::Matrix<ComplexOrRealType>& a_,
                          const PsimagLite::Matrix<ComplexOrRealType>& b_,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin_,
                          SizeType offsetY,
                          typename PsimagLite::Vector<ComplexOrRealType>::Type& xout_,
                          SizeType offsetX)
{
	const bool is_complex = PsimagLite::IsComplexNumber<ComplexOrRealType>::True;
	const int nrow_A = a_.n_row();
	const int ncol_A = a_.n_col();
	const int nrow_B = b_.n_row();
	const int ncol_B = b_.n_col();

	const int isTransA = (transA == 'T') || (transA == 't');
	const int isTransB = (transB == 'T') || (transB == 't');
	const int isConjTransA = (transA == 'C') || (transA == 'c');
	const int isConjTransB = (transB == 'C') || (transB == 'c');

	const int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	const int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;
	const int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	const int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

	const int nrow_X = nrow_2;
	const int ncol_X = nrow_1;
	const int nrow_Y = ncol_2;
	const int ncol_Y = ncol_1;

	PsimagLite::MatrixNonOwned<ComplexOrRealType> xout(nrow_X, ncol_X, xout_, offsetX);
	PsimagLite::MatrixNonOwned<const ComplexOrRealType> yin(nrow_Y, ncol_Y, yin_, offsetY);

	assert((imethod == 1) ||
	       (imethod == 2) ||
	       (imethod == 3));

/*
 *   -------------------------------------------------------------
 *   A and B in dense matrix format
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




	if (imethod == 1) {

	/*
	 *  --------------------------------------------
	 *  BY(iby,jby) = op(B(ib,jb))*Y(iy,jy)
	 *
	 *  X(ix,jx) += BY(iby,jby ) * transpose(op(A(ia,ja)))
	 *  --------------------------------------------
	 */
		const int nrow_BY = nrow_X;
		const int ncol_BY = ncol_Y;
		PsimagLite::Matrix<ComplexOrRealType> by_(nrow_BY, ncol_BY );
		PsimagLite::MatrixNonOwned<ComplexOrRealType> byRef(by_);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> byConstRef(by_);

	/*
	 * ---------------
	 * setup BY(ib,ja)
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
	 * BY(iby,jby)  = op(B(ib,jb))*Y(iy,jy)
	 * ------------------------------
	 */
			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;
			den_matmul_pre(  trans,
			                 nrow_B,
			                 ncol_B,
			                 b_,

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
	 * X(ix,jx) += BY(iby,jby) * transpose(op(A(ia,ja)))
	 * -------------------------------------------
	 */
			const char trans = (isTransA || isConjTransA) ? 'N' : 'T';
			if (is_complex && isConjTransA) {
				// --------------------------------------------
				// transpose( conj( transpose(A) ) ) is conj(A)
				// perform  conj operation
				// --------------------------------------------
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
		const int nrow_YAt = nrow_Y;
		const int ncol_YAt = ncol_X;
		PsimagLite::Matrix<ComplexOrRealType> yat_(nrow_YAt, ncol_YAt);
		PsimagLite::MatrixNonOwned<ComplexOrRealType> yatRef(yat_);
		PsimagLite::MatrixNonOwned<const ComplexOrRealType> yatConstRef(yat_);

	/*
	 * ----------------
	 * setup YAt(jb,ia)
	 * ----------------
	 */




		{
			int iy = 0;
			int jy = 0;

			// not needed, FIXME
			for(jy=0; jy < ncol_YAt; jy++) {
				for(iy=0; iy < nrow_YAt; iy++) {
					yat_(iy,jy) = 0;
				};
			};
		}


		{
			/*
	  * ---------------------
	  * YAt(jb,ia) = Y(jb,ja) * tranpose(op(A(ia,ja)))
	  * ---------------------
	  */
			const char trans = (isTransA || isConjTransA) ? 'N' : 'T';
			if (is_complex && isConjTransA) {
				// --------------------------------------------
				// transpose( conj( transpose(A) ) ) is conj(A)
				// perform in-place conj operation
				// --------------------------------------------
                                PsimagLite::Matrix<ComplexOrRealType> a_conj(nrow_A,ncol_A);

				for(int ja=0; ja < ncol_A; ja++) {
				for(int ia=0; ia < nrow_A; ia++) {
					a_conj(ia,ja) = PsimagLite::conj( a_(ia,ja) );
				};
				};

			        den_matmul_post( trans,
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



			        den_matmul_post( trans,
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
	 * X(ib,ia) += op(B(ib,jb)) * YAt(jb,ia)
	 * ------------
	 */

			// const char trans = (isTransB) ? 'T' : 'N';
			const char trans = transB;
			den_matmul_pre( trans,
			                nrow_B,
			                ncol_B,
			                b_,

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
	*
	* C = kron(transpose(A),B)
	* C([ib,ja], [jb,ia]) = At(ja,ia) * B(ib,jb)
	* X([ib,ja]) = B(ib,jb) * Y(jb,ia) * transpose(At(ja,ia))
	* X([ib,ja]) = B(ib,jb) * Y(jb,ia) * A(ia,ja)
	*            = (A(ia,ja)*B(ib,jb)) * Y(jb,ia)
	*
	* C = kron(A, transpose(B))
	* C([jb,ia],[ib,ja]) = A(ia,ja) * Bt(jb,ib)
	* X(jb,ia) = (A(ia,ja) * Bt(jb,ib)) * Y(ib,ja)
	* X(jb,ia) = Bt(jb,ib) * Y(ib,ja) * transpose(A(ia,ja))
	*          = Bt(jb,ib) * Y(ib,ja) * At(ja,ia)
	*          = B(ib,jb)  * Y(ib,ja) * A(ia,ja)
	*          = (A(ia,ja)*B(ib,jb)) * Y(ib,ja)
	*
	*
	* C = kron( transpose(A), transpose(B))
	* C([jb,ja], [ib,ia] ) = At(ja,ia) * Bt(jb,ib)
	* X(jb,ja) = ( At(ja,ia) * Bt(jb,ib) ) * Y(ib,ia)
	*          = Bt(jb,ib) * Y(ib,ia) * transpose(At(ja,ia))
	*          = B(ib,jb) * Y(ib,ia) * A(ia,ja)
	* ---------------------------------------------
	*/

		int ia = 0;
		int ja = 0;
		int ib = 0;
		int jb = 0;

		for(ia=0; ia < nrow_A; ia++) {
			for(ja=0; ja < ncol_A; ja++) {
				for(ib=0; ib < nrow_B; ib++) {
					for(jb=0; jb < ncol_B; jb++) {
						ComplexOrRealType aij = a_(ia,ja);
						if (is_complex && isConjTransA) {
							aij = PsimagLite::conj( aij );
						};

						ComplexOrRealType bij = b_(ib,jb);
						if (is_complex && isConjTransB) {
							bij = PsimagLite::conj( bij );
						};

						ComplexOrRealType cij = aij * bij;

						int ix = (isTransB || isConjTransB) ? jb : ib;
						int jx = (isTransA || isConjTransA) ? ja : ia;
						int iy = (isTransB || isConjTransB) ? ib : jb;
						int jy = (isTransA || isConjTransA) ? ia : ja;

						ComplexOrRealType yij = yin(iy,jy);
						xout(ix,jx) +=  (cij * yij);
					};
				};
			};
		};

	};

}

template<typename ComplexOrRealType>
void den_kron_mult(const char transA,
                   const char transB,
                   const PsimagLite::Matrix<ComplexOrRealType>& a_,
                   const PsimagLite::Matrix<ComplexOrRealType>& b_,
                   const typename PsimagLite::Vector<ComplexOrRealType>::Type& yin,
                   SizeType offsetY,
                   typename PsimagLite::Vector<ComplexOrRealType>::Type& xout,
                   SizeType offsetX,
                   const typename PsimagLite::Real<ComplexOrRealType>::Type denseFlopDiscount)
{
/*
 *   -------------------------------------------------------------
 *   A and B in dense matrix format
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
	const int nrow_A = a_.n_row();
	const int ncol_A = a_.n_col();
	const int nrow_B = b_.n_row();
	const int ncol_B = b_.n_col();
	int nnz_A = nrow_A * ncol_A;
	int nnz_B = nrow_B * ncol_B;

	ComplexOrRealType kron_nnz = 0;
	ComplexOrRealType kron_flops = 0;
	int imethod = 1;

	const int isTransA = (transA == 'T') || (transA == 't');
	const int isTransB = (transB == 'T') || (transB == 't');
	const int isConjTransA = (transA == 'C') || (transA == 'c');
	const int isConjTransB = (transB == 'C') || (transB == 'c');

	int nrow_1 = (isTransA || isConjTransA) ? ncol_A : nrow_A;
	int ncol_1 = (isTransA || isConjTransA) ? nrow_A : ncol_A;


	int nrow_2 = (isTransB || isConjTransB) ? ncol_B : nrow_B;
	int ncol_2 = (isTransB || isConjTransB) ? nrow_B : ncol_B;

	estimate_kron_cost( nrow_1,ncol_1,nnz_A, nrow_2,ncol_2,nnz_B,
	                    &kron_nnz, &kron_flops, &imethod, denseFlopDiscount);



	den_kron_mult_method(
	            imethod,

	            transA,
	            transB,

	            a_,

	            b_,

	            yin,
	            offsetY,
	            xout,
	            offsetX);



}
#undef BY
#undef YAt
#undef X
#undef Y
#undef A
#undef B
#undef X2
