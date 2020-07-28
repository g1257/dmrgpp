#ifndef DMFTSOLVER_H
#define DMFTSOLVER_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"

namespace Dmft {

template<typename ComplexOrRealType>
class DmftSolver {

public:

	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;
	typedef FunctionOfFrequencyType::RealType RealType;
	typedef Dmft::Dispersion DispersionType;

	DmftSolver(RealType fictiousBeta, SizeType nMatsubaras, const DispersionType& dispersion)
	    : sigma_(fictiousBeta, nMatsubaras),
	      latticeG_(fictiousBeta, nMatsubaras),
	      dispersion_(dispersion)
	{}

	void selfConsistencyLoop()
	{
		computeLatticeGf();
	}

private:

	void computeLatticeGf()
	{
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		SizeType totalKvalues = dispersion_.size();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType wn = sigma_.omega(i);
			const ComplexOrRealType value = sigma_(i);
			for (SizeType j = 0; j < totalKvalues; ++j ) {
				RealType ek = dispersion_(j);
				sum += 1.0/(wn - ek + mu_ - value);
			}

			latticeG_(i) = sum/totalKvalues;
		}
	}

	FunctionOfFrequencyType sigma_;
	FunctionOfFrequencyType latticeG_;
	const DispersionType& dispersion_;
};
}
#endif // DMFTSOLVER_H
