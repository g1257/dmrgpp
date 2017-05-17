#include "util.h"

template<typename ComplexOrRealType>
void den_gen_matrix(const int nrow_A,
                    const int ncol_A,
                    const typename PsimagLite::Real<ComplexOrRealType>::Type threshold,
                    PsimagLite::Matrix<ComplexOrRealType>& a_)
{
	/*
 * -------------------------------
 * generate a random matix in (0,1)
 * accept only if   aij < threshold
 * full matrix if threshold > 1
 * sparse matrix if threshold << 1
 * -------------------------------
 */

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	const ComplexOrRealType dzero = 0;

	int ia = 0;
	int ja = 0;

	for(ja=0; ja < ncol_A; ja++) {
		for(ia=0; ia < nrow_A; ia++) {
			RealType drand = rand()/static_cast<RealType>(RAND_MAX);
			RealType aij   = rand()/static_cast<RealType>(RAND_MAX);

			int is_accept = (drand <= threshold);

			a_(ia,ja) = (is_accept) ? aij : dzero;

		};
	};
}
#undef A
