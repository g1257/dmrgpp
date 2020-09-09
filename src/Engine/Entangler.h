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
		for (SizeType i = 0; i < n; ++i) {
			if (seen[i]) continue;
			seen[i] = true;
			WordType ket = basis[i];
			WordType bket = bar(ket);
			if (ket == bket) { // in class C
				f(i, i) = 1;
				continue;
			}

			auto it = std::find(basis.begin(), basis.end(), bket);
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
