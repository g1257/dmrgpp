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
	typedef typename FunctionOfFrequencyType::MatsubarasType MatsubarasType;
	typedef typename MatsubarasType::VectorRealType VectorRealType;
	typedef Dmft::Dispersion<ComplexOrRealType> DispersionType;
	typedef Fit<ComplexOrRealType> FitType;
	typedef typename FitType::MinParamsType MinParamsType;
	typedef ParamsDmftSolver<ComplexOrRealType, InputNgType> ParamsDmftSolverType;
	typedef ImpuritySolver<ParamsDmftSolverType> ImpuritySolverType;
	typedef typename ImpuritySolverType::ApplicationType ApplicationType;

	DmftSolver(const ParamsDmftSolverType& params, const ApplicationType& app)
	    : params_(params),
	      sigma_(params.ficticiousBeta, params.nMatsubaras),
	      latticeG_(params.ficticiousBeta, params.nMatsubaras),
	      gammaG_(params.ficticiousBeta, params.nMatsubaras),
	      dispersion_(params.numberOfKpoints),
	      mu_(params.mu),
	      fit_(params.nBath, params.minParams),
	      impuritySolver_(params, app)
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

		if (error < params_.dmftError) {
			std::cout<<"Converged after "<<iter<<" iterations; error="<<error<<"\n";
			return; // <--- EARLY EXIT HERE
		}

		std::cout<<"I did "<<iter<<" iterations; but error="<<error;
		std::cout<<" is greater than the tolerance="<<params_.dmftError;
		std::cout<<" that was requested\n";
	}

	void print(std::ostream& os) const
	{
		os<<"Sigma\n";
		os<<sigma_;
		printBathParams(os);
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
			const ComplexOrRealType diff = oldValue - newValue;
			sum += PsimagLite::real(diff*PsimagLite::conj(diff));
			sigma_(i) = newValue;
		}

		return sum;
	}

	void printBathParams(std::ostream& os) const
	{
		os<<"bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath\n";
		os<<"bathParams[nBath-...] ==> energies on each bath site\n";
		os<<fit_.result();
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
