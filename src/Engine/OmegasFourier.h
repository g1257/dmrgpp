#ifndef OMEGASFOURIER_H
#define OMEGASFOURIER_H

namespace Dmrg {

template<typename RealType>
class OmegasFourier {

public:

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	void fourier(VectorRealType& qValues,
	             const VectorRealType& values1,
	             const VectorRealType& values2)
	{
		err("unimplemented fourier\n");
	}


	void writeFourier(VectorRealType& array, const VectorRealType& qValues)
	{
		err("unimplemented writeFourier\n");
	}
};
}
#endif // OMEGASFOURIER_H
