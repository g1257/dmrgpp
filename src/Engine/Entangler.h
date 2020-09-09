#ifndef ENTANGLER_H
#define ENTANGLER_H
#include "Vector.h"
#include "Matrix.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename HilbertBasisType, typename ComplexOrRealType>
class Entangler {

public:

	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename HilbertBasisType::value_type WordType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	/* F^+ | alpha, \bar{alpha} > = |\bar{alpha}, \alpha> if alpha in A or C.
	                              = 0 otherwise

	We partition the states of one site in three disjoint classes A, B, C, as follows.
	We put all states of one-site in a bucket. We take a state alpha out of the bucket,
	and also \bar{alpha} out of the bucket. If \bar{alpha} = \alpha we put alpha in C.
	If \bar{alpha} \neq  \alpha, we put alpha in A and \bar{alpha} in B.
	We proceed like this until there are no more states left in the bucket.
	*/

	template<typename BarType>
	static void setGammaMatrix(MatrixType& f,
	                           const HilbertBasisType& basis,
	                           const BarType& bar)
	{
		const SizeType n = basis.size();
		VectorBoolType seen(n, false);
		SizeType nroot = sqrt(n);
		if (nroot*nroot != n)
			err("Basis is not a perfect square\n");
		nroot = ProgramGlobals::logBase2(nroot);
		WordType mask = (1 << nroot) - 1;

		for (SizeType i = 0; i < n; ++i) {
			if (seen[i]) continue;
			seen[i] = true;
			WordType ket = basis[i];

			WordType physKet = ket & mask;
			ket >>= nroot;
			WordType ancKet = ket & mask;

			WordType physBar = bar(physKet);
			if (physBar != ancKet) continue; // this is a "bad" state

			if (physBar == physKet) { // in class C
				f(i, i) = 1;
				continue;
			}

			// this is not needed because it's taken care of
			// by the seen variable
			//if (physBar < physKet) continue; // in class B

			// in class A
			WordType barFull = (physKet << nroot);
			barFull |= physBar;

			auto it = std::find(basis.begin(), basis.end(), barFull);
			assert(it != std::end(basis));
			SizeType j = it - basis.begin();
			assert(j < seen.size());
			if (seen[j]) continue;
			seen[j] = true;
			f(i, j) = 1;
		}
	}
};
}
#endif // ENTANGLER_H
