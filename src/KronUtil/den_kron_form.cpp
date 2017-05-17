#include "util.h"

template<typename ComplexOrRealType>
void den_kron_form( const int nrow_A, 
                    const int ncol_A,
                    const PsimagLite::Matrix<ComplexOrRealType>& a_,
                    const int nrow_B,
                    const int ncol_B,
                    const PsimagLite::Matrix<ComplexOrRealType>& b_,
                    PsimagLite::Matrix<ComplexOrRealType>& c_ )
{
	/*
 * ---------------------------------------
 * form C = kron(A,B),  where
 * nrow_C = nrow_A * nrow_B
 * ncol_C = ncol_A * ncol_B
 * C([ib,ia], [jb,ja]) = A(ia,ja)*B(ib,jb)
 * ---------------------------------------
 */
	int ia = 0;
	int ja = 0;
	int ib = 0;
	int jb = 0;


	for(ja=0; ja < ncol_A; ja++) {
		for(jb=0; jb < ncol_B; jb++) {
			for(ia=0; ia < nrow_A; ia++) {
				for(ib=0; ib < nrow_B; ib++) {
					int ic = ib + ia*nrow_B;
					int jc = jb + ja*ncol_B;

					c_(ic,jc) = a_(ia,ja) * b_(ib,jb);
				};
			};
		};
	};
}
