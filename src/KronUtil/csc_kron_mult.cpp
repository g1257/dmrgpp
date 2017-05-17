#include "util.h"

template<typename ComplexOrRealType>
void csc_kron_mult_method(const int imethod,
                          const int nrow_A,
                          const int ncol_A,
                          const PsimagLite::Vector<int>::Type& acolptr,
                          const PsimagLite::Vector<int>::Type& arow,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                          const int nrow_B,
                          const int ncol_B,
                          const PsimagLite::Vector<int>::Type& bcolptr,
                          const PsimagLite::Vector<int>::Type& brow,
                          const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,
                          const PsimagLite::Matrix<ComplexOrRealType>& yin,
                          PsimagLite::Matrix<ComplexOrRealType>& xout)
{
	const int nrow_X = nrow_B;
	const int ncol_X = nrow_A;

	const int nrow_Y = ncol_B;
	const int ncol_Y = ncol_A;

	int nnz_A = csc_nnz( ncol_A, acolptr );
	int nnz_B = csc_nnz( ncol_B, bcolptr );
	int has_work = (nnz_A >= 1) && (nnz_B >= 1);

	assert((imethod == 1) ||
	       (imethod == 2) ||
	       (imethod == 3));

	if (!has_work) {
		return;
	};
	/*
 *   -------------------------------------------------------------
 *   A and B are in compressed sparse COLUMN format
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
	 *  BY(ib,ja) = B(ib,jb)*Y(jb,ja)
	 *
	 *  X(ib,ia) += BY(ib,ja ) * transpose(A(ia,ja))
	 *  --------------------------------------------
	 */



		/*
	 * ---------------
	 * setup BY(ib,ja)
	 * ---------------
	 */
		const int nrow_BY = nrow_B;
		const int ncol_BY = ncol_A;
		PsimagLite::Matrix<ComplexOrRealType> by_(nrow_BY, ncol_BY);

		{
			int iby = 0;
			int jby = 0;

			// not needed, FIXME
			for( jby=0; jby < ncol_BY; jby++) {
				for( iby=0; iby < nrow_BY; iby++) {
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
			const char trans = 'N';
			csc_matmul_pre(  trans,
			                 nrow_B,
			                 ncol_B,
			                 bcolptr,
			                 brow,
			                 bval,

			                 nrow_Y,
			                 ncol_Y,
			                 yin,

			                 nrow_BY,
			                 ncol_BY,
			                 by_);

		}

		{
			/*
	 * -------------------------------------------
	 * X(ib,ia) += BY(ib,ja) * transpose(A(ia,ja))
	 * -------------------------------------------
	 */
			const char trans = 'T';

			csc_matmul_post(
			            trans,
			            nrow_A,
			            ncol_A,
			            acolptr,
			            arow,
			            aval,

			            nrow_BY,
			            ncol_BY,
			            by_,

			            nrow_X,
			            ncol_X,
			            xout);

		}
	}
	else if (imethod == 2) {
		/*
	 * ---------------------
	 * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja))
	 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
	 * ---------------------
	 */



		/*
	 * ----------------
	 * setup YAt(jb,ia)
	 * ----------------
	 */

		PsimagLite::Matrix<ComplexOrRealType> yat_(ncol_B, nrow_A);

		int nrow_YAt = ncol_B;
		int ncol_YAt = nrow_A;


		{
			int jb = 0;
			int ia = 0;

			//not needed, FIXME
			for(ia=0; ia < ncol_A; ia++) {
				for(jb=0; jb < ncol_B; jb++) {
					yat_(jb,ia) = 0;
				};
			};
		}


		{
			/*
	  * ---------------------
	  * YAt(jb,ia) = Y(jb,ja) * tranpose(A(ia,ja)
	  * ---------------------
	  */
			const char transa = 'T';



			csc_matmul_post( transa,
			                 nrow_A,
			                 ncol_A,
			                 acolptr,
			                 arow,
			                 aval,

			                 nrow_Y,
			                 ncol_Y,
			                 yin,

			                 nrow_YAt,
			                 ncol_YAt,
			                 yat_);
		}




		{
			/*
	 * ------------
	 * X(ib,ia) += B(ib,jb) * YAt(jb,ia)
	 * ------------
	 */

			const char trans = 'N';
			csc_matmul_pre( trans,
			                nrow_B,
			                ncol_B,
			                bcolptr,
			                brow,
			                bval,

			                nrow_YAt,
			                ncol_YAt,
			                yat_,

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

		int ja = 0;
		for(ja=0; ja < ncol_A; ja++) {
			int istarta = acolptr[ja];
			int ienda = acolptr[ja+1]-1;
			int ka = 0;
			for(ka=istarta; ka <= ienda; ka++) {
				int ia = arow[ka];
				ComplexOrRealType aij = aval[ka];
				int jb = 0;
				for(jb=0; jb < ncol_B; jb++) {
					int istartb = bcolptr[jb];
					int iendb = bcolptr[jb+1]-1;
					int kb = 0;
					for(kb=istartb; kb <= iendb; kb++) {
						int ib = brow[kb];
						ComplexOrRealType bij = brow[kb];
						ComplexOrRealType cij = aij * bij;

						xout(ib,ia) +=  (cij * yin(jb,ja));
					};
				};
			};
		};

	};

}

template<typename ComplexOrRealType>
void csc_kron_mult( const int nrow_A,
                    const int ncol_A,
                    const PsimagLite::Vector<int>::Type& acolptr,
                    const PsimagLite::Vector<int>::Type& arow,
                    const typename PsimagLite::Vector<ComplexOrRealType>::Type& aval,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Vector<int>::Type& bcolptr,
                    const PsimagLite::Vector<int>::Type& brow,
                    const typename PsimagLite::Vector<ComplexOrRealType>::Type& bval,
                    const PsimagLite::Matrix<ComplexOrRealType>& yin,
                    PsimagLite::Matrix<ComplexOrRealType>& xout )

{
	/*
 *   -------------------------------------------------------------
 *   A and B are in compressed sparse COLUMN format
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

	int nnz_A = csc_nnz( ncol_A, acolptr );
	int nnz_B = csc_nnz( ncol_B, bcolptr );
	int has_work = (nnz_A >= 1) && (nnz_B >= 1);

	ComplexOrRealType kron_nnz = 0;
	ComplexOrRealType kron_flops = 0;
	int imethod = 1;

	if (!has_work) {
		return;
	};

	estimate_kron_cost( nrow_A,ncol_A,nnz_A, nrow_B,ncol_B,nnz_B,
	                    &kron_nnz, &kron_flops, &imethod );



	csc_kron_mult_method(
	            imethod,
	            nrow_A,
	            ncol_A,
	            acolptr,
	            arow,
	            aval,

	            nrow_B,
	            ncol_B,
	            bcolptr,
	            brow,
	            bval,

	            yin,
	            xout );
}
#undef BY
#undef YAt
#undef X
#undef Y
