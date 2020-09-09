#ifndef ENTANGLER_H
#define ENTANGLER_H
#include "Vector.h"
#include "Matrix.h"

namespace Dmrg {

template<typename HilbertBasisType, typename ComplexOrRealType>
class Entangler {

public:

	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename HilbertBasisType::value_type WordType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	template<typename BarType>
	static void setGammaMatrix(MatrixType& f,
	                           const HilbertBasisType& basis,
	                           const BarType& bar)
	{
		const SizeType n = basis.size();
		VectorBoolType seen(n);
		const SizeType nroot = sqrt(n);
		if (nroot*nroot != n)
			err("Basis is not a perfect square\n");
		WordType mask = (1 << nroot) - 1;

		for (SizeType i = 0; i < n; ++i) {
			if (seen[i]) continue;
			seen[i] = true;
			WordType ket = basis[i];

			WordType physKet = ket & mask;
			ket >>= nroot;
			WordType ancKet = ket & mask;

			WordType physBar = bar(physKet);
			if (physKet != ancKet) continue; // this is a "bad" state

			if (physBar == physKet) { // in class C
				f(i, i) = 1;
				continue;
			}

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
