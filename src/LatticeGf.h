#ifndef LATTICEGF_H
#define LATTICEGF_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"
#include "PsimagLite.h"

namespace Dmft {

template<typename ComplexOrRealType>
class LatticeGf {

public:

	class DensityOfStates {

	public:

		DensityOfStates(PsimagLite::String option)
		{
			err("DensityOfStates not yet supported\n");
		}
	};

	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;
	typedef Dispersion<ComplexOrRealType> DispersionType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef DensityOfStates DensityOfStatesType;

	LatticeGf(const FunctionOfFrequencyType& sigma,
	          RealType mu,
	          PsimagLite::String option)
	    : dispersion_(nullptr),
	      dos_(nullptr),
	      sigma_(sigma),
	      mu_(mu),
	      latticeG_(sigma.fictitiousBeta(), sigma.totalMatsubaras()),
	      gammaG_(sigma.fictitiousBeta(), sigma.totalMatsubaras())
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, option, ",");
		if (tokens.size() == 0)
			err("LatticeGf: INTERNAL ERROR: No option?!\n");

		PsimagLite::String option0 = tokens[0];
		PsimagLite::String option1 = (tokens.size() >= 2) ? tokens[1] : "";
		PsimagLite::String option2 = (tokens.size() >= 3) ? tokens[2] : "";

		if (tokens.size() > 3)
			printCoutCerr("LatticeGf: Ignoring extra options after " + tokens[2] + "\n");

		if (option0 == "momentum") {
			const SizeType kpoints = PsimagLite::atoi(option2);
			dispersion_ = new DispersionType(option1, kpoints);
		} else if (option0 == "energy") {
			dos_ = new DensityOfStatesType(option1);
			if (option2 != "")
				printCoutCerr("LatticeGf: Ignoring extra option" + option2 + "\n");
		} else {
			err("LatticeGf: Unknow option " + option0 +
			    ". Only momentum or energy allowed\n");
		}
	}

	~LatticeGf()
	{
		delete dispersion_;
		dispersion_ = nullptr;

		delete dos_;
		dos_ = nullptr;
	}

	const FunctionOfFrequencyType& operator()() const { return latticeG_; }

	const FunctionOfFrequencyType& gammaG() const { return gammaG_; }

	void update()
	{
		if (dispersion_)
			updateMomentum();
		else
			err("DOS:: Unimplemented\n");
	}

private:

	void updateMomentum() {

		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		SizeType totalKvalues = dispersion_->size();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType value = sigma_(i);
			ComplexOrRealType sum = 0.0;
			for (SizeType j = 0; j < totalKvalues; ++j ) {
				RealType ek = dispersion_->operator ()(j);
				sum += 1.0/(iwn -ek + mu_ - value);
			}

			latticeG_(i) = sum/static_cast<RealType>(totalKvalues);
			gammaG_(i) = iwn - 1.0/latticeG_(i) - value;
		}
	}

	static void printCoutCerr(PsimagLite::String str)
	{
		std::cout<<str;
		std::cerr<<str;
	}

	LatticeGf(const LatticeGf&) = delete;

	LatticeGf& operator=(const LatticeGf&) = delete;

	DispersionType* dispersion_;
	DensityOfStatesType* dos_;
	const FunctionOfFrequencyType& sigma_;
	RealType mu_;
	FunctionOfFrequencyType latticeG_;
	FunctionOfFrequencyType gammaG_;
};

}
#endif // LATTICEGF_H
