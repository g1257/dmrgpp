#ifndef DMFTSOLVER_H
#define DMFTSOLVER_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"
#include "Fit.h"
#include "ParamsDmftSolver.h"
#include "ImpuritySolver.h"

namespace Dmft {

template<typename ComplexOrRealType, typename InputNgType>
class DmftSolver {

public:

	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;
	typedef typename FunctionOfFrequencyType::RealType RealType;
	typedef typename FunctionOfFrequencyType::VectorRealType VectorRealType;
	typedef Dmft::Dispersion<ComplexOrRealType> DispersionType;
	typedef Fit<ComplexOrRealType> FitType;
	typedef typename FitType::MinParamsType MinParamsType;
	typedef ParamsDmftSolver<ComplexOrRealType, InputNgType> ParamsDmftSolverType;
	typedef ImpuritySolver<ParamsDmftSolverType> ImpuritySolverType;

	DmftSolver(const ParamsDmftSolverType& params)
	    : params_(params),
	      sigma_(params.ficticiousBeta, params.nMatsubaras),
	      latticeG_(params.ficticiousBeta, params.nMatsubaras),
	      gammaG_(params.ficticiousBeta, params.nMatsubaras),
	      dispersion_(params.numberOfKpoints),
	      mu_(params.mu),
	      fit_(params.nBath, params.minParams),
	      impuritySolver_(params)
	{}

	// DMFT Self consistency loop; see Steve Johnston's notes
	void selfConsistencyLoop()
	{
		SizeType iter = 0;
		RealType error = 0;

		const SizeType totalMatsubaras = sigma_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) sigma_(i) = 1.0/(i+1.0);

		for (; iter < params_.dmftIter; ++iter) {

			std::cout<<"SelfConsistLoop iter= "<<iter<<"\n";
			computeLatticeGf();

			fit_.fit(gammaG_);

			impuritySolver_.solve(fit_.result());

			error = computeNewSelfEnergy(fit_.result());

			std::cout<<"SelfConsistLoop error="<<error<<"\n";
			if (error < params_.dmftError)
				break;
		}

		std::cout<<"Converged after "<<iter<<" iterations; error="<<error<<"\n";
	}

private:

	void computeLatticeGf()
	{
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		SizeType totalKvalues = dispersion_.size();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType value = sigma_(i);
			ComplexOrRealType sum = 0.0;
			for (SizeType j = 0; j < totalKvalues; ++j ) {
				RealType ek = dispersion_(j);
				sum += 1.0/(iwn -ek + mu_ - value);
			}

			latticeG_(i) = sum/static_cast<RealType>(totalKvalues);
			gammaG_(i) = iwn - 1.0/latticeG_(i) - value;
		}
	}

	RealType computeNewSelfEnergy(const VectorRealType& bathParams)
	{
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		RealType sum = 0;
		typename FitType::AndersonFunctionType andersonFunction(params_.nBath, gammaG_);
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType oldValue = sigma_(i);
			const ComplexOrRealType newValue = iwn - andersonFunction.anderson(bathParams, i) -
			        1.0/impuritySolver_.gimp(i);
			sum += PsimagLite::real(FitType::AndersonFunctionType::squareOf(oldValue - newValue));
			sigma_(i) = newValue;
		}

		return sum;
	}

	const ParamsDmftSolverType& params_;
	FunctionOfFrequencyType sigma_;
	FunctionOfFrequencyType latticeG_;
	FunctionOfFrequencyType gammaG_;
	DispersionType dispersion_;
	RealType mu_;
	FitType fit_;
	ImpuritySolverType impuritySolver_;
};
}
#endif // DMFTSOLVER_H
