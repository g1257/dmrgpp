#ifndef IMAGINARYUNITORFAIL_HH
#define IMAGINARYUNITORFAIL_HH
#include "Complex.h"
#include "PsimagLite.h"

namespace Dmrg
{

template <typename T>
class ImginaryUnitOrFail
{
public:

	static T value()
	{
		err("This features can only be run when usecomplex is in SolverOptions\n");
		return 0.;
	}
};

template <typename T>
class ImginaryUnitOrFail<std::complex<T>>
{
public:

	static std::complex<T> value()
	{
		return std::complex<T>(0., 1.);
	}
};
}
#endif // IMAGINARYUNITORFAIL_HH
