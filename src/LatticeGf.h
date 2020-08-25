#ifndef LATTICEGF_H
#define LATTICEGF_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"
#include "PsimagLite.h"
#include "Integrator.h"

namespace Dmft {

template<typename ComplexOrRealType>
class LatticeGf {

	class DensityOfStates {

	public:

		typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

		DensityOfStates(PsimagLite::String option, RealType wOverTwo)
		    : wOverTwo_(wOverTwo)
		{
			if (option != "semicircular")
				err("DensityOfStates " + option +
				    " not yet supported; only semicircular supported.\n");

		}

		RealType lowerBound() const { return -wOverTwo_; }

		RealType upperBound() const { return wOverTwo_; }

		RealType operator()(RealType e) const
		{
			const RealType wOverTwoSquared = wOverTwo_*wOverTwo_;
			return 2.0*sqrt(wOverTwoSquared - e*e)/(wOverTwoSquared*M_PI);
		}

	private:

		RealType wOverTwo_;
	};

	class Integrand {

		struct Params {

			Params(DensityOfStates* dos_, ComplexOrRealType iwnMinusSigma_)
			    : dos(dos_), iwnMinusSigma(iwnMinusSigma_)
			{}

			Params(const Params&) = delete;

			Params& operator=(const Params&) = delete;

			DensityOfStates* dos;
			ComplexOrRealType iwnMinusSigma;
		};

	public:

		typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

		Integrand(DensityOfStates* dos, ComplexOrRealType iwnMinusSigma)
		    : p_(dos, iwnMinusSigma)
		{}

		static RealType function(RealType x, void* vp)
		{
			Params* p = static_cast<Params*>(vp);
			return p->dos->operator()(x)/(p->iwnMinusSigma - x);
		}

		void update(ComplexOrRealType iwnMinusSigma)
		{
			p_.iwnMinusSigma_ = iwnMinusSigma;
		}

		//Params& params() { return p_; }

	private:

		Params p_;
	};

public:

	typedef FunctionOfFrequency<ComplexOrRealType> FunctionOfFrequencyType;
	typedef Dispersion<ComplexOrRealType> DispersionType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

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
			const RealType W = PsimagLite::atof(option2);
			dos_ = new DensityOfStates(option1, 0.5*W);
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

	void updateEnergy()
	{
		Integrand integrand(dos_, 0.0);
		PsimagLite::Integrator<Integrand> integrator(integrand);
		typename PsimagLite::Vector<RealType>::Type pts(2,0);
		pts[0] = dos_->lowerBound();
		pts[1] = dos_->upperBound();
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType value = sigma_(i);
			integrand.update(iwn - value);
			latticeG_(i) = integrator(pts);
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
	DensityOfStates* dos_;
	const FunctionOfFrequencyType& sigma_;
	RealType mu_;
	FunctionOfFrequencyType latticeG_;
	FunctionOfFrequencyType gammaG_;
};

}
#endif // LATTICEGF_H
