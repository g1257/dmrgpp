#ifndef IMPURITYSOLVER_BASE_H
#define IMPURITYSOLVER_BASE_H

#include "CincuentaInputCheck.h"
#include "FreqEnum.h"
#include "InputNg.h"
#include "Matsubaras.h"
#include "ModelParams.h"
#include "PsimagLite.h"
#include "RealFrequencyRange.hh"
#include "Vector.h"

namespace Dmft {

template <typename ComplexOrRealType> class ImpuritySolverBase {

public:

	using RealType               = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType            = std::complex<RealType>;
	using VectorRealType         = typename PsimagLite::Vector<RealType>::Type;
	using VectorComplexType      = typename PsimagLite::Vector<ComplexType>::Type;
	using ApplicationType        = PsimagLite::PsiApp;
	using ModelParamsType        = ModelParams<RealType>;
	using InputNgType            = PsimagLite::InputNg<CincuentaInputCheck>;
	using InputNgReadableType    = InputNgType::Readable;
	using MatsubarasType         = PsimagLite::Matsubaras<RealType>;
	using RealFrequencyRangeType = PsimagLite::RealFrequencyRange<RealType>;

	// ctor
	ImpuritySolverBase(RealType                                            ficticiousBeta,
	                   SizeType                                            nMatsubaras,
	                   PsimagLite::InputNg<CincuentaInputCheck>::Readable& io)
	    : matsubaras_(ficticiousBeta, nMatsubaras, 0.) // last argument is real part
	    , real_range_(io)
	{ }

	// Virtuals BEGIN
	// dtor
	virtual ~ImpuritySolverBase() { }

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	virtual void
	solve(const VectorRealType& bathParams, PsimagLite::FreqEnum freq_enum, SizeType iter)
	    = 0;

	virtual const VectorComplexType& gimp() const = 0;

	virtual PsimagLite::FreqEnum freqEnum() const = 0;

	virtual const MatsubarasType& matsubaras() const { return matsubaras_; }

	virtual const RealFrequencyRangeType& realFreqRange() const { return real_range_; }
	// Virtuals END

	// public static START
	static ComplexOrRealType density(const VectorComplexType& g)
	{
		const SizeType    n   = g.size();
		ComplexOrRealType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += g[i];

		return sum;
	}

	static void scale(std::vector<ComplexOrRealType>& g)
	{
		const ComplexOrRealType factor = 0.25;
		for (auto& val : g) {
			val *= factor;
		}
	}
	// public static END

protected:

	enum class GsOrOmegaEnum
	{
		GS,
		OMEGA
	};

	static std::string createGsInput(const ModelParamsType& model_params,
	                                 InputNgReadableType&   io)
	{
		std::string        data  = commonInputString(model_params, io, GsOrOmegaEnum::GS);
		PsimagLite::String data2 = addBathParams(data, model_params);
		return data2;
	}

	static std::string commonInputString(const ModelParamsType& model_params,
	                                     InputNgReadableType&   io,
	                                     GsOrOmegaEnum          gs_or_omega)
	{
		RealType U = 0;
		io.readline(U, "HubbardU=");
		std::string hubbardU_vector = buildHubbardU(U, model_params.numberOfSites());

		std::string additional_solver_options;
		try {
			io.readline(additional_solver_options, "SolverOptions=");
		} catch (std::exception&) { }

		if (gs_or_omega == GsOrOmegaEnum::OMEGA) {
			additional_solver_options
			    += ",CorrectionVectorTargeting,restart,minimizeDisk";
		}

		std::string root;
		io.readline(root, "RootOutputname=");
		std::string gs_output    = root + "gs";
		std::string omega_output = root + "omega";
		std::string output = (gs_or_omega == GsOrOmegaEnum::GS) ? gs_output : omega_output;

		SizeType infinite_loops = 0;
		io.readline(infinite_loops, "InfiniteLoopKeptStates=");

		std::string label = "FiniteLoops";
		if (gs_or_omega == GsOrOmegaEnum::GS) {
			label += "Gs=";
		} else {
			label += "Omega=";
		}

		std::string finite_loops;
		io.readline(finite_loops, label);

		SizeType nup = 0;
		io.readline(nup, "TargetElectronsUp=");

		SizeType ndown = 0;
		io.readline(ndown, "TargetElectronsDown=");

		std::string s
		    = "##Ainur1.0\n\nTotalNumberOfSites=" + ttos(model_params.numberOfSites())
		    + ";\nNumberOfTerms=1;\nDegreesOfFreedom=1;\nGeometryKind=star;"
		      "\nGeometryOptions=none;\nhubbardU="
		    + hubbardU_vector
		    + ";\nModel=HubbardOneBand;\nSolverOptions=twositedmrg,geometryallinsystem,"
		      "hd5dontprint";
		if (!additional_solver_options.empty()) {
			s += "," + additional_solver_options;
		}

		s += ";\nVersion=templateForDMFT;\nOutputFile=" + output
		    + ";\nInfiniteLoopKeptStates=" + ttos(infinite_loops) + ";\n";
		s += "FiniteLoops=" + finite_loops + ";\nTargetElectronsUp=" + ttos(nup)
		    + ";\nTargetElectronsDown=" + ttos(ndown) + ";\n";

		return s;
	}

	static std::string addBathParams(const std::string&     data,
	                                 const ModelParamsType& model_params)
	{
		const PsimagLite::String connectors = vectorToString(model_params.hoppings(), ",");
		const PsimagLite::String label      = "dir0:Connectors=[" + connectors + "];\n";
		const PsimagLite::String potentialV
		    = vectorToString(model_params.potentialV(), ",");
		const PsimagLite::String label2
		    = "potentialV=[" + potentialV + "," + potentialV + "];\n";

		return data + label + label2;
	}

private:

	static std::string buildHubbardU(const RealType& U, SizeType n)
	{
		std::string s = "[" + ttos(U);
		for (SizeType i = 1; i < n; ++i) {
			s += ", 0.";
		}
		s += "]";
		return s;
	}

	static PsimagLite::String vectorToString(const VectorRealType& v,
	                                         const std::string&    separator)
	{
		PsimagLite::String buffer;
		SizeType           n = v.size();
		for (SizeType i = 0; i < n; ++i) {
			buffer += ttos(v[i]);
			if (i + 1 < n) {
				buffer += ",";
			}
		}

		return buffer;
	}

	MatsubarasType         matsubaras_;
	RealFrequencyRangeType real_range_;
};
}
#endif // IMPURITYSOLVER_BASE_H
