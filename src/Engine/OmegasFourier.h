#ifndef OMEGASFOURIER_H
#define OMEGASFOURIER_H

namespace Dmrg {

template<typename RealType, typename ReadableType>
class OmegasFourier {

public:

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	OmegasFourier(bool skipFourier) : skipFourier_(skipFourier)
	{}

	void configure(ReadableType& io)
	{
		if (skipFourier_) return;
	}

	void fourier(const VectorRealType& values1,
	             const VectorRealType& values2)
	{
		if (skipFourier_) return;
		//use qvalues_;
		err("unimplemented fourier\n");
	}

	void writeFourier()
	{
		if (skipFourier_) return;
		//use qvalues_;
		//VectorRealType array;
		err("unimplemented writeFourier\n");
	}

	void printGnuplot()
	{
		if (skipFourier_) return;
		err("unimplemented writeFourier\n");
	}

private:

	bool skipFourier_;
	VectorRealType qValues_;
};
}
#endif // OMEGASFOURIER_H
