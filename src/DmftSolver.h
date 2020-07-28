#ifndef DMFTSOLVER_H
#define DMFTSOLVER_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"
#include "Fit.h"
#include "ParamsDmftSolver.h"

namespace Dmft {

template<typename ComplexOrRealType, typename InputNgType>
class DmftSolver {

public:

	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;
	typedef typename FunctionOfFrequencyType::RealType RealType;
	typedef Dmft::Dispersion<ComplexOrRealType> DispersionType;
	typedef Fit<ComplexOrRealType> FitType;
	typedef typename FitType::MinParamsType MinParamsType;
	typedef ParamsDmftSolver<ComplexOrRealType, InputNgType> ParamsDmftSolverType;

	DmftSolver(const ParamsDmftSolverType& params)
	    : params_(params),
	      sigma_(params.ficticiousBeta, params.nMatsubaras),
	      latticeG_(params.ficticiousBeta, params.nMatsubaras),
	      gammaG_(params.ficticiousBeta, params.nMatsubaras),
	      dispersion_(params.numberOfKpoints),
	      mu_(params.mu),
	      fit_(params.nBath, params.minParams)
	{}

	void selfConsistencyLoop()
	{
		computeLatticeGf();

		fit_.fit();
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

	const ParamsDmftSolverType& params_;
	FunctionOfFrequencyType sigma_;
	FunctionOfFrequencyType latticeG_;
	FunctionOfFrequencyType gammaG_;
	DispersionType dispersion_;
	RealType mu_;
	FitType fit_;
};
}
#endif // DMFTSOLVER_H
