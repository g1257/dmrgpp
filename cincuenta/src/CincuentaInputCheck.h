/*! \file CincuentaInputCheck.h
 *
 *  CincuentaInputChecking functions
 */
#ifndef CINCUENTA_INPUT_CHECK_H
#define CINCUENTA_INPUT_CHECK_H
#include "InputCheck.h"
#include <PsimagLite/Geometry/Geometry.h>
#include <PsimagLite/Options.h>
#include <stdexcept>
#include <vector>
// #include "ProgramGlobals.h"

namespace Dmft {

class CincuentaInputCheck {

	using OptionsReadableType = PsimagLite::Options::Readable;
	using VectorStringType    = PsimagLite::Vector<PsimagLite::String>::Type;

public:

	CincuentaInputCheck()
	    : optsReadable_(0)
	{
		// knownLabels_.push_back("TotalNumberOfSites");
	}

	~CincuentaInputCheck()
	{
		if (optsReadable_ != 0)
			delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		// PsimagLite::String str = PsimagLite::Geometry<int,int,ProgramGlobals>::import();
		PsimagLite::String str = Dmrg::InputCheck().import();

		str += "integer FicticiousBeta;\n";
		str += "real ChemicalPotential;\n";
		str += "integer Matsubaras;\n";
		str += "integer NumberOfBathPoints;\n";
		str += "integer DmftNumberOfIterations;\n";
		str += "real DmftTolerance;\n";
		str += "real MinParamsDelta;\n";
		str += "real MinParamsDelta2;\n";
		str += "real MinParamsTolerance;\n";
		str += "integer MinParamsMaxIter;\n";
		str += "integer MinParamsVerbose;\n";
		str += "string LatticeGf;\n";
		str += "string ImpuritySolver;\n";
		str += "real InitBathRa;\n";
		str += "real InitBathRb;\n";
		str += "string FitMethod;\n";
		str += "string RootOutputname;\n";
		str += "string FitOptions;\n";

		// Non-equilibrium DMFT (interaction quench)
		str += "real HubbardUFinal;\n";
		str += "real TmaxNeq;\n";
		str += "integer NtNeq;\n";
		str += "integer NeqDmftIter;\n";
		// Keep the removed labels declared solely so InputNg can reach the
		// migration rejection in checkSimpleLabel(), rather than issuing its
		// generic undeclared-label parse error first.
		str += "string NeqSolver;\n";
		str += "real NeqDmftTolerance;\n";
		// GBEK two-bath scheme (Gramsch, Balzer, Eckstein, Kollar PRB 88, 235106).
		str += "integer NeqBathRank;\n"; // L: rank of Cholesky second bath
		str += "real BandwidthFinal;\n"; // W_f for hopping quench; default 0 = no quench
		str += "string QuenchShape;\n"; // "step" (default), "cosine", "tanh"
		str += "real QuenchDuration;\n"; // ramp duration t_q in real time; 0 = step quench
		str += "string NeqOutputPrefix;\n"; // prefix for output Green's function files
		str += "integer NeqAtomicLimit;\n"; // 1 = start neq run from the atomic limit:
		                                    // no first bath (Delta^- = 0 exactly),
		                                    // per GBEK Sec. VI setup

		return str;
	}

	bool check(const PsimagLite::String&                           label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType                                            line) const
	{
		return false;
	}

	bool
	check(const PsimagLite::String& label, const PsimagLite::String& vec, SizeType line) const
	{
		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label, SizeType line) const
	{
		// InputNg calls this while reading every label, before any typed
		// readline() call. Reject removed labels here so they cannot be merely
		// warned about and then silently ignored by the GBEK-only dispatcher.
		if (label == "NeqSolver")
			return error1(
			    "NeqSolver= was removed: GBEK is now the only non-equilibrium "
			    "solver. Remove NeqSolver= and configure GBEK with "
			    "NeqBathRank>0 and NeqDmftIter>0",
			    line);
		if (label == "NeqDmftTolerance")
			return error1("NeqDmftTolerance= was removed: GBEK uses the fixed "
			              "corrector count NeqDmftIter, not a convergence tolerance. "
			              "Remove NeqDmftTolerance= and set NeqDmftIter>0",
			              line);

		for (SizeType i = 0; i < knownLabels_.size(); ++i)
			if (knownLabels_[i] == label)
				return true;
		PsimagLite::String msg("WARNING: Unknown label " + label + "\n");
		std::cout << msg;
		std::cerr << msg;
		return false;
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr << "USAGE is " << name << "\n";
	}

	// Ainur validates typed labels before InputNg invokes checkSimpleLabel().
	// Every executable using this input check must call this immediately after
	// constructing its Readable input, so typed legacy labels receive the same
	// useful migration error as legacy input syntax.
	template <typename ReadableType> static void rejectRemovedLabels(ReadableType& io)
	{
		PsimagLite::String solver;
		if (hasLabel(io, solver, "NeqSolver="))
			throw PsimagLite::RuntimeError(
			    "NeqSolver= was removed: GBEK is now the only non-equilibrium "
			    "solver; remove NeqSolver= and configure GBEK with NeqBathRank>0 "
			    "and NeqDmftIter>0\n");

		double tolerance = 0;
		if (hasLabel(io, tolerance, "NeqDmftTolerance="))
			throw PsimagLite::RuntimeError(
			    "NeqDmftTolerance= was removed: GBEK uses the fixed corrector "
			    "count NeqDmftIter, not a convergence tolerance; remove "
			    "NeqDmftTolerance= and set NeqDmftIter>0\n");
	}

private:

	template <typename ReadableType, typename T>
	static bool hasLabel(ReadableType& io, T& value, const PsimagLite::String& label)
	{
		try {
			io.readline(value, label);
		} catch (std::exception&) {
			return false;
		}
		return true;
	}

	bool checkForVector(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() == 0)
			return false;
		SizeType n = atoi(vec[0].c_str());
		return (vec.size() == n + 1);
	}

	bool checkForMatrix(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() < 2)
			return false;
		SizeType row = atoi(vec[0].c_str());
		SizeType col = atoi(vec[1].c_str());
		SizeType n   = row * col;
		return (vec.size() == n + 2);
	}

	bool error1(const PsimagLite::String& message, SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	OptionsReadableType* optsReadable_;
	VectorStringType     allowedFileOptions_;
	VectorStringType     knownLabels_;
}; // class CincuentaInputCheck
} // namespace Dmft

/*@}*/
#endif
