#ifndef DMFTSOLVER_H
#define DMFTSOLVER_H
#include "Fit.h"
#include "FunctionOfFrequency.h"
#include "ImpuritySolverBase.h"
#include "ImpuritySolverDmrg.h"
#include "ImpuritySolverExactDiag.h"
#include "LatticeGf.h"
#include "ParamsDmftSolver.h"
#include <PsimagLite/InputNg.h>
#include <fstream>
#include <iomanip>
#include <limits>

namespace Dmft {

template <typename ComplexOrRealType> class DmftSolver {

public:

	using FunctionOfFrequencyType     = FunctionOfFrequency<ComplexOrRealType>;
	using RealType                    = typename FunctionOfFrequencyType::RealType;
	using MatsubarasType              = typename FunctionOfFrequencyType::MatsubarasType;
	using RealFrequencyRangeType      = PsimagLite::RealFrequencyRange<RealType>;
	using VectorRealType              = typename MatsubarasType::VectorRealType;
	using FitType                     = Fit<ComplexOrRealType>;
	using MinParamsType               = typename FitType::MinParamsType;
	using ParamsDmftSolverType        = ParamsDmftSolver<ComplexOrRealType>;
	using ImpuritySolverType          = ImpuritySolverBase<ComplexOrRealType>;
	using ImpuritySolverExactDiagType = ImpuritySolverExactDiag<ComplexOrRealType>;
	using ImpuritySolverDmrgType      = ImpuritySolverDmrg<ComplexOrRealType>;
	using LatticeGfType               = LatticeGf<ComplexOrRealType>;
	using ApplicationType             = typename ImpuritySolverType::ApplicationType;
	using AndersonFunctionType        = typename FitType::AndersonFunctionType;
	using VectorComplexType           = typename ImpuritySolverType::VectorComplexType;
	using InputNgType                 = PsimagLite::InputNg<CincuentaInputCheck>;

	// Equilibrium information needed by a future neq initialization. It owns
	// the Matsubara data because the subsequent real-frequency solve replaces
	// impuritySolver_->gimp().
	struct EquilibriumInitialData {
		VectorRealType    bathParameters;
		VectorComplexType gimpMatsubara;
		VectorRealType    matsubaraFrequencies;

		void writeMatsubara(const std::string& filename) const
		{
			if (matsubaraFrequencies.size() != gimpMatsubara.size())
				err("EquilibriumInitialData: frequency/value size mismatch\n");

			std::ofstream output(filename);
			if (!output)
				err("EquilibriumInitialData: cannot open " + filename + "\n");

			output << std::setprecision(std::numeric_limits<RealType>::max_digits10);
			for (SizeType i = 0; i < matsubaraFrequencies.size(); ++i) {
				const ComplexOrRealType value = gimpMatsubara[i];
				output << matsubaraFrequencies[i] << " " << PsimagLite::real(value)
				       << " " << PsimagLite::imag(value) << "\n";
			}

			if (!output)
				err("EquilibriumInitialData: failed while writing " + filename
				    + "\n");
		}
	};

	DmftSolver(const ParamsDmftSolverType&          params,
	           const typename FitType::InitResults& initResults,
	           const ApplicationType&               app,
	           InputNgType::Readable&               io)
	    : params_(params)
	    , sigma_(params.ficticiousBeta, params.nMatsubaras)
	    , latticeG_(sigma_, params.mu, params.latticeGf)
	    , fit_(params.nBath, params.minParams, initResults)
	    , impuritySolver_(nullptr)
	    , io_(io)
	{
		if (params.impuritySolver == "dmrg")
			impuritySolver_ = new ImpuritySolverDmrgType(params, app, io_);
		else if (params.impuritySolver == "exactdiag")
			impuritySolver_ = new ImpuritySolverExactDiagType(params, app, io_);
		else
			err("Unknown impurity solver " + params.impuritySolver + "\n");

		// FitOptions="particleholesymmetric" makes AndersonFunction fit mirror-image
		// bath pairs (epsilon, -epsilon) about ChemicalPotential=; that bath is only
		// the correct particle-hole-symmetric one if ChemicalPotential also sits at
		// the impurity's actual particle-hole-symmetric point, which the solved
		// impurity model fixes at U/2 (see ModelParams.h). A ChemicalPotential=
		// away from U/2 would fit a symmetric-looking bath to the wrong target.
		if (FitType::computeOptions(params.fit_options)
		    == FitType::Options::PARTICLE_HOLE_SYMM) {
			RealType U = 0;
			io.readline(U, "HubbardU=");
			const RealType muExpected = RealType(0.5) * U;
			if (std::abs(params.mu - muExpected) > 1e-10)
				err("DmftSolver: FitOptions=\"particleholesymmetric\" requires "
				    "ChemicalPotential=HubbardU/2 ("
				    + ttos(muExpected)
				    + "); got ChemicalPotential=" + ttos(params.mu) + "\n");
		}
	}

	~DmftSolver()
	{
		delete impuritySolver_;
		impuritySolver_ = nullptr;
	}

	// DMFT Self consistency loop; see Steve Johnston's notes
	void selfConsistencyLoop()
	{
		SizeType iter  = 0;
		RealType error = 0;

		// const SizeType totalMatsubaras = sigma_.totalMatsubaras();
		// for (SizeType i = 0; i < totalMatsubaras; ++i) sigma_(i) = 1.0/(i+1.0);

		typename FitType::Options fit_options
		    = FitType::computeOptions(params_.fit_options);

		for (; iter < params_.dmftIter; ++iter) {

			std::cout << "SelfConsistLoop iter= " << iter << "\n";

			latticeG_.update();

			fit_.fit(latticeG_.g0(), params_.mu, fit_options);

			impuritySolver_->solve(
			    fit_.result(), PsimagLite::FreqEnum::MATSUBARA, iter);

			// For now this is a dependable way to
			// make sure the we get normalized results
			// forom the impuritySolver. This could be
			// removed if all impurity solvers gave G with
			// the expected weighting.
			impuritySolver_->enforceSpectralSumRule();

			this->logDebug();

			error = computeNewSelfEnergy(fit_.result());

			std::cout << "SelfConsistLoop error=" << error << "\n";

			if (error < params_.dmftError)
				break;
		}

		snapshotEquilibriumInitialData();

		impuritySolver_->solve(fit_.result(), PsimagLite::FreqEnum::REAL, 0);
		this->logDebug();

		if (error < params_.dmftError) {
			std::cout << "Converged after " << iter << " iterations; error=" << error
			          << "\n";
			return; // <--- EARLY EXIT HERE
		}

		std::cout << "I did " << iter << " iterations; but error=" << error;
		std::cout << " is greater than the tolerance=" << params_.dmftError;
		std::cout << " that was requested\n";
	}

	const VectorRealType& bathResult() const { return fit_.result(); }

	const EquilibriumInitialData& equilibriumInitialData() const
	{
		return equilibriumInitialData_;
	}

	void print(std::ostream& os) const
	{
		os << "Sigma\n";
		os << sigma_;

		printBathParams(os);

		writeGimpForDebugOnly(os);

		os << "LatticeG\n";
		os << latticeG_();

		os << "G0\n";
		os << latticeG_.g0();

		FunctionOfFrequencyType siteEx(sigma_.fictitiousBeta(),
		                               sigma_.totalMatsubaras() / 2);

		AndersonFunctionType anderson_function(params_.nBath, params_.mu);
		computeSiteExcludedG(siteEx, anderson_function);
		os << "SiteExcludedG\n";
		os << siteEx;

		printAndersonFunction(os, anderson_function);
		printClusterG0(os, anderson_function);
	}

private:

	void snapshotEquilibriumInitialData()
	{
		assert(impuritySolver_->freqEnum() == PsimagLite::FreqEnum::MATSUBARA);
		equilibriumInitialData_.bathParameters = fit_.result();
		equilibriumInitialData_.gimpMatsubara  = impuritySolver_->gimp();

		const MatsubarasType& matsubaras      = impuritySolver_->matsubaras();
		const SizeType        totalMatsubaras = matsubaras.total();
		assert(equilibriumInitialData_.gimpMatsubara.size() == totalMatsubaras);
		equilibriumInitialData_.matsubaraFrequencies.resize(totalMatsubaras);
		for (SizeType i = 0; i < totalMatsubaras; ++i)
			equilibriumInitialData_.matsubaraFrequencies[i] = matsubaras.omega(i);
	}

	void printAndersonFunction(std::ostream& os, const AndersonFunctionType& af) const
	{
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		os << "AndersonFunction\n";
		os << totalMatsubaras << "\n";
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const RealType          wn = sigma_.omega(i);
			const ComplexOrRealType val
			    = af.anderson(fit_.result(), ComplexOrRealType(0, wn));
			os << wn << " " << val << "\n";
		}
	}

	void printClusterG0(std::ostream& os, const AndersonFunctionType& af) const
	{
		SizeType totalMatsubaras = sigma_.totalMatsubaras();
		os << "Gcluster0\n";
		os << totalMatsubaras << "\n";
		ClusterG0<ComplexOrRealType> cluster_g0(af);
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			RealType          wn = sigma_.omega(i);
			ComplexOrRealType val
			    = cluster_g0.g0cluster(fit_.result(), ComplexOrRealType(0, wn));
			os << wn << " " << val << "\n";
		}
	}

	void writeGimpForDebugOnly(const std::string& file_out) const
	{
		std::ofstream fout(file_out);
		if (!fout || !fout.good())
			err(std::string("Could not write to") + file_out + "\n");

		writeGimpForDebugOnly(fout);

		fout.close();
	}

	void writeGimpForDebugOnly(std::ostream& os) const
	{
		if (impuritySolver_->freqEnum() == PsimagLite::FreqEnum::MATSUBARA) {
			writeGimpDebugOnly(os, impuritySolver_->matsubaras());
		} else {
			writeGimpDebugOnly(os, impuritySolver_->realFreqRange());
		}
	}

	template <typename SomeFreqRangeType>
	void writeGimpDebugOnly(std::ostream& os, const SomeFreqRangeType& freq_range) const
	{
		const VectorComplexType& gimp = impuritySolver_->gimp();
		const SizeType           n    = gimp.size();

		for (SizeType i = 0; i < n; ++i) {
			const ComplexOrRealType value = gimp[i];
			const RealType          omega = freq_range.omega(i);
			os << omega << " " << PsimagLite::real(value) << " "
			   << PsimagLite::imag(value) << "\n";
		}
	}

	void logDebug() const
	{
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);

		if (mpiRank != 0) {
			return;
		}

		// temporary to fix legacy name
		std::remove("gimp_exact.txt");

		const VectorComplexType& gimp      = impuritySolver_->gimp();
		PsimagLite::FreqEnum     freq_enum = impuritySolver_->freqEnum();

		ComplexOrRealType d = ImpuritySolverType::density(gimp);
		std::cerr << "Sum of Gimp=" << d << "\n";

		std::string root = "gimp_" + params_.impuritySolver;

		if (freq_enum == PsimagLite::FreqEnum::MATSUBARA) {
			writeGimpForDebugOnly(root + ".txt");
			writeLatticeGForDebug("latticeG_" + params_.impuritySolver + ".txt");
		} else {
			writeGimpForDebugOnly(root + "_real.txt");
		}
	}

	// Write the local lattice GF to a file.
	// Format: one line per Matsubara frequency  "omega Re Im"  (no header count).
	// For U=0 (no interactions, Sigma=0), this must equal gimp to within bath fitting error.
	void writeLatticeGForDebug(const std::string& filename) const
	{
		std::ofstream fout(filename);
		if (!fout || !fout.good())
			err(std::string("Could not write to ") + filename + "\n");
		const FunctionOfFrequencyType& g = latticeG_();
		const SizeType                 n = g.totalMatsubaras();
		for (SizeType i = 0; i < n; ++i)
			fout << g.omega(i) << " " << PsimagLite::real(g(i)) << " "
			     << PsimagLite::imag(g(i)) << "\n";
	}

	RealType computeNewSelfEnergy(const VectorRealType& bathParams)
	{
		SizeType                               totalMatsubaras = sigma_.totalMatsubaras();
		RealType                               sum             = 0;
		typename FitType::AndersonFunctionType andersonFunction(params_.nBath, params_.mu);

		const VectorComplexType& gimp = impuritySolver_->gimp();
		assert(gimp.size() == totalMatsubaras);
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			const ComplexOrRealType iwn      = ComplexOrRealType(0.0, sigma_.omega(i));
			const ComplexOrRealType oldValue = sigma_(i);
			const ComplexOrRealType newValue = iwn + params_.mu
			    - andersonFunction.anderson(bathParams, iwn) - 1.0 / gimp[i];
			const ComplexOrRealType diff = oldValue - newValue;
			sum += PsimagLite::real(diff * PsimagLite::conj(diff));
			sigma_(i) = newValue;
		}

		return sum / totalMatsubaras;
	}

	void computeSiteExcludedG(FunctionOfFrequencyType&    siteEx,
	                          const AndersonFunctionType& anderson_function) const
	{
		const SizeType totalMatsubaras = siteEx.totalMatsubaras();
		for (SizeType i = 0; i < totalMatsubaras; ++i) {
			ComplexOrRealType iwn = ComplexOrRealType(0.0, siteEx.omega(i));
			ComplexOrRealType sumOverAlpha
			    = anderson_function.anderson(fit_.result(), iwn);
			ComplexOrRealType value = iwn - sumOverAlpha;
			siteEx(i)               = 1.0 / value;
		}
	}

	void printBathParams(std::ostream& os) const
	{
		os << "bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath\n";
		os << "bathParams[nBath-...] ==> energies on each bath site\n";
		os << fit_.result();
	}

	const ParamsDmftSolverType& params_;
	FunctionOfFrequencyType     sigma_;
	LatticeGfType               latticeG_;
	FitType                     fit_;
	ImpuritySolverType*         impuritySolver_;
	EquilibriumInitialData      equilibriumInitialData_;
	InputNgType::Readable&      io_;
};
}
#endif // DMFTSOLVER_H
