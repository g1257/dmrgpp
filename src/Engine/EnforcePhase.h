#ifndef ENFORCEPHASE_H
#define ENFORCEPHASE_H
#include "Matrix.h"

namespace Dmrg {

template<typename ComplexOrRealType>
class EnforcePhase {

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

public:

	static void enforcePhase(PsimagLite::Matrix<ComplexOrRealType>& a)
	{
		SizeType cols = a.cols();
		for (SizeType i = 0; i < cols; ++i) {
			enforcePhase(a, i);
		}
	}

private:

	static void enforcePhase(PsimagLite::Matrix<ComplexOrRealType>& a, SizeType col)
	{
		RealType sign1 = 0;
		SizeType rows = a.rows();
		for (SizeType j = 0; j < rows; ++j) {
			RealType val = PsimagLite::norm(a(j, col));
			if (val < 1e-6) continue;
			sign1 = (val > 0) ? 1 : -1;
			break;
		}

		assert(sign1 != 0);
		// get a consistent phase
		for (SizeType j = 0; j < rows; ++j)
			a(j, col) *= sign1;
	}
};
}
#endif // ENFORCEPHASE_H
