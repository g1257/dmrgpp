#ifndef LATTICEGF_H
#define LATTICEGF_H
#include "FunctionOfFrequency.h"
#include "Dispersion.h"
#include "DensityOfStates.h"
#include "PsimagLite.h"
#include "Integrator.h"

namespace Dmft {

template<typename ComplexOrRealType>
class LatticeGf {

	typedef DensityOfStates<ComplexOrRealType> DensityOfStatesType;

	template<int RealOrImg>
	class Integrand {

		struct Params {

			Params(DensityOfStatesType* dos_, ComplexOrRealType iwnMinusSigma_)
			    : dos(dos_), iwnMinusSigma(iwnMinusSigma_)
			{}

			Params(const Params&) = delete;

			Params& operator=(const Params&) = delete;

			DensityOfStatesType* dos;
			ComplexOrRealType iwnMinusSigma;
		};

	public:

		typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

		Integrand(DensityOfStatesType* dos, ComplexOrRealType iwnMinusSigma)
		    : p_(dos, iwnMinusSigma)
		{}

		static RealType function(RealType x, void* vp)
		{
			Params* p = static_cast<Params*>(vp);
			ComplexOrRealType result = p->dos->operator()(x)/(p->iwnMinusSigma - x);
			return (RealOrImg == 0) ? PsimagLite::real(result) : PsimagLite::imag(result);
		}

		void update(ComplexOrRealType iwnMinusSigma)
		{
			p_.iwnMinusSigma = iwnMinusSigma;
		}

		Params& params() { return p_; }

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
	      latticeG_(sigma.fictitiousBeta(), sigma.totalMatsubaras()/2),
	      gammaG_(sigma.fictitiousBeta(), sigma.totalMatsubaras()/2)
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
			std::cout<<"LatticeGf: Using momentum with " + option1 + " and " +
			           option2 + " k points.\n";
		} else if (option0 == "energy") {
			const RealType W = PsimagLite::atof(option2);
			dos_ = new DensityOfStatesType(option1, 0.5*W);
			std::cout<<"LatticeGf: Using energy with " + option1 + " and W = " +
			           option2 + "\n";
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
			updateEnergy();
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
		Integrand<0> integrand0(dos_, 0.0); // real part
		PsimagLite::Integrator<Integrand<0> > integrator0(integrand0);

		Integrand<1> integrand1(dos_, 0.0); // imag part
		PsimagLite::Integrator<Integrand<1> > integrator1(integrand1);

		typename PsimagLite::Vector<RealType>::Type pts(2,0);
		pts[0] = dos_->lowerBound();
		pts[1] = dos_->upperBound();
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType value = sigma_(i);
			integrand0.update(iwn - value);
			integrand1.update(iwn - value);
			latticeG_(i) = ComplexOrRealType(integrator0(pts), integrator1(pts));
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
