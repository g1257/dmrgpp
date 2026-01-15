#ifndef PEIERLSHOPPING_HH
#define PEIERLSHOPPING_HH
#include "Complex.h"
#include <utility>
#include <vector>

namespace Dmrg {

template <typename SuperGeometryType, bool> class BuildPierls {
public:

	using ComplexOrRealType = typename SuperGeometryType::ComplexOrRealType;
	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType    = std::vector<RealType>;

	static auto lambda(const VectorRealType& A)
	{
		err("Cannot run Peierls Hubbard with reals; say usecomplex in SolverOptions\n");
		return [A](ComplexOrRealType&, RealType, SizeType) { };
	}
};

template <typename SuperGeometryType> class BuildPierls<SuperGeometryType, true> {
public:

	using ComplexOrRealType = typename SuperGeometryType::ComplexOrRealType;
	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using VectorRealType    = std::vector<RealType>;

	static auto lambda(const VectorRealType& A)
	{
		return [&A_const
		        = std::as_const(A)](std::complex<RealType>& hop, RealType t, SizeType site)
		{
			if (A_const.size() == 0)
				return;
			assert(site < A_const.size());
			RealType arg = A_const[site] * t;
			hop *= std::complex<RealType>(cos(arg), sin(arg));
		};
	}
};

}
#endif // PEIERLSHOPPING_HH
