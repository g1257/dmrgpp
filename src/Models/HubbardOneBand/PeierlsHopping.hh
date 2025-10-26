#ifndef PEIERLSHOPPING_HH
#define PEIERLSHOPPING_HH
#include "Complex.h"

namespace Dmrg {

template<typename ComplexOrRealType>
class BuildPierls
{
public:
	using RealType = typename PsimagLite::Real<ComplexOrRealType>::Type;

	static auto lambda(const RealType& A)
	{
		err("Cannot run Peierls Hubbard with reals; say usecomplex in SolverOptions\n");
                return [A](ComplexOrRealType&, RealType) {};
	}
};

template<typename RealType>
class BuildPierls<std::complex<RealType>>
{
public:

	static auto lambda(const RealType& A)
	{
		return [A](std::complex<RealType>& hop, RealType t)
		{
			RealType arg = A * t;
			hop *= std::complex<RealType>(cos(arg), sin(arg));
		};
	}
};

}
#endif // PEIERLSHOPPING_HH
