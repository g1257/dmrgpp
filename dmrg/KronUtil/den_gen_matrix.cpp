#include "util.h"

#include <random>

template <typename ComplexOrRealType>
void den_gen_matrix(const int                                                 nrow_A,
                    const int                                                 ncol_A,
                    const typename PsimagLite::Real<ComplexOrRealType>::Type& threshold,
                    PsimagLite::Matrix<ComplexOrRealType>&                    a_)
{
	/*
	 * -------------------------------
	 * generate a random matix in (0,1)
	 * accept only if   aij < threshold
	 * full matrix if threshold > 1
	 * sparse matrix if threshold << 1
	 * -------------------------------
	 */

	using RealType                = typename PsimagLite::Real<ComplexOrRealType>::Type;
	const ComplexOrRealType dzero = 0;

	std::random_device rd;
	std::mt19937       rng(rd());

	std::uniform_real_distribution<RealType> value(-1.0, 1.0);
	std::uniform_real_distribution<RealType> keep(0.0, 1.0);

	for (int ja = 0; ja < ncol_A; ja++) {
		for (int ia = 0; ia < nrow_A; ia++) {
			if constexpr (PsimagLite::IsComplexNumber<ComplexOrRealType>::True)
				a_(ia, ja) = (keep(rng) <= threshold)
				    ? ComplexOrRealType(value(rng), value(rng))
				    : dzero;
			else
				a_(ia, ja) = (keep(rng) <= threshold) ? value(rng) : dzero;
		}
	}
}
