#ifndef DMFTSOLVER_H
#define DMFTSOLVER_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"

namespace Dmft {

template<typename ComplexOrRealType>
class DmftSolver {

public:

	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;
	typedef typename FunctionOfFrequencyType::RealType RealType;
	typedef Dmft::Dispersion<ComplexOrRealType> DispersionType;

	DmftSolver(RealType fictiousBeta,
	           SizeType nMatsubaras,
	           const DispersionType& dispersion,
	           RealType mu)
	    : sigma_(fictiousBeta, nMatsubaras),
	      latticeG_(fictiousBeta, nMatsubaras),
	      gammaG_(fictiousBeta, nMatsubaras),
	      dispersion_(dispersion),
	      mu_(mu)
	{}

	void selfConsistencyLoop()
	{
		computeLatticeGf();

		fitBathParams();

	}

private:

	void computeLatticeGf()
	{
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		SizeType totalKvalues = dispersion_.size();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType wn = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType value = sigma_(i);
			ComplexOrRealType sum = 0.0;
			for (SizeType j = 0; j < totalKvalues; ++j ) {
				RealType ek = dispersion_(j);
				sum += 1.0/(wn -ek + mu_ - value);
			}

			latticeG_(i) = sum/static_cast<RealType>(totalKvalues);
			gammaG_(i) = wn - 1.0/latticeG_(i) - value;
		}
	}

	void fitBathParams()
	{

	}

	FunctionOfFrequencyType sigma_;
	FunctionOfFrequencyType latticeG_;
	FunctionOfFrequencyType gammaG_;
	const DispersionType& dispersion_;
	RealType mu_;
};
}
#endif // DMFTSOLVER_H
