/*
Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[LanczosPlusPlus++, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include "Geometry/Geometry.h"
#include "Options.h"
#include <stdexcept>
#include <vector>

namespace LanczosPlusPlus {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;

public:

	InputCheck()
	    : optsReadable_(0)
	{ }

	~InputCheck()
	{
		if (optsReadable_ != 0)
			delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		PsimagLite::String str = PsimagLite::Geometry<int, int, LanczosGlobals>::import();

		str += "vector hubbardU;\n";
		str += "vector potentialV;\n";
		str += "string! Model;\n";
		str += "string! SolverOptions;\n";
		str += "string! Version;\n";
		str += "integer! InfiniteLoopKeptStates;\n";
		str += "string! OutputFile;\n";
		str += "matrix.integer FiniteLoops;\n";
		str += "integer RepeatFiniteLoopsFrom;\n";
		str += "integer RepeatFiniteLoopsTimes;\n";
		str += "integer TargetElectronsUp;\n";
		str += "integer TargetElectronsDown;\n";
		str += "integer TargetElectronsTotal;\n";
		str += "real GsWeight;\n";
		str += "real TSPTau;\n";
		str += "integer TSPTimeSteps;\n";
		str += "integer TSPAdvanceEach;\n";
		str += "string TSPAlgorithm;\n";
		str += "vector.integer TSPSites;\n";
		str += "vector.integer TSPLoops;\n";
		str += "string TSPProductOrSum;\n";
		str += "string TSPOperator;\n";
		str += "string OperatorExpression;\n";
		str += "integer Threads = 1;\n";
		str += "integer Orbitals = 1;\n";
		str += "string FeAsMode;\n";
		str += "integer TargetSpinTimesTwo;\n";
		str += "integer UseSu2Symmetry;\n";
		str += "integer Pvectors;\n";
		str += "string TruncationTolerance;\n";
		str += "integer HeisenbergTwiceS;\n";
		str += "integer TargetSzPlusConst;\n";
		str += "integer SpinTwiceS;\n";
		str += "integer OrbitalTwiceS;\n";
		str += "real LambdaOne;\n";
		str += "real LambdaTwo;\n";
		str += "real CorrectionA;\n";
		str += "string RestartFilename;\n";
		str += "real LanczosEps;\n";
		str += "real TridiagonalEps;\n";
		str += "integer DynamicDmrgType;\n";
		str += "real CorrectionVectorFreqType;\n";
		str += "real CorrectionVectorEta;\n";
		str += "string CorrectionVectorAlgorithm;\n";
		str += "real CorrectionVectorOmega;\n";
		str += "string Intent;\n";
		str += "integer OpOnSiteThreshold;\n";
		str += "integer FirstRitz;\n";
		str += "integer CVnForFraction;\n";
		str += "real AnisotropyD;\n";
		str += "string FindSymmetrySector;\n";
		str += "string AddOnSiteHamiltonian;\n";
		str += "vector MagneticFieldX;\n";
		str += "vector MagneticFieldZ;\n";

		return str;
	}

	bool check(const PsimagLite::String&                           label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType                                            line) const
	{
		if (label == "JMVALUES") {
			if (vec.size() != 2)
				return error1("JMVALUES", line);
			return true;
		} else if (label == "RAW_MATRIX" || label == "SpinOrbit") {
			SizeType row = atoi(vec[0].c_str());
			SizeType col = atoi(vec[1].c_str());
			SizeType n   = row * col;
			if (vec.size() != n + 2)
				return error1("RAW_MATRIX", line);
			return true;
		} else if (label == "Connectors") {
			return true;
		} else if (label == "MagneticField") {
			return true;
		} else if (label == "FiniteLoops") {
			SizeType n = atoi(vec[0].c_str());
			if (vec.size() != 3 * n + 1)
				return error1("FiniteLoops", line);
			return true;
		}
		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label, SizeType line) const
	{
		// FIXME: needs implementation
		return true;
	}

	void check(const PsimagLite::String& label, const PsimagLite::String& val, SizeType)
	{
		if (label != "SolverOptions")
			return;
		PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

		/* PSIDOC LanczosSolverOptions
		\begin{itemize}
		\item[none] Use this as a placeholder. ``none'' does not disable other options.
		\item[InternalProductStored] Stored the sparse matrix in memory before diagonalizing
		it.
		\item[InternalProductOnTheFly] Compute the sparse matrix on-the-fly while
		diagonalizing it.
		\item[printmatrix] Print the Hamiltonian matrix.
		\item[dumpmatrix] Use exact diagonalization instead of Lanczos diagonalization,
		and output all information to obtain the full spectrum.
		\item [setAffinities] TBW
		\end{itemize}
		*/
		registerOpts.push_back("none");
		registerOpts.push_back("InternalProductStored");
		registerOpts.push_back("InternalProductOnTheFly");
		registerOpts.push_back("printmatrix");
		registerOpts.push_back("dumpmatrix");
		registerOpts.push_back("setAffinities");

		PsimagLite::Options::Writeable optWriteable(
		    registerOpts, PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new OptionsReadableType(optWriteable, val);
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	void usage(const char* progName)
	{
		std::cerr << "Usage: " << progName << " [-g -c] -f filename\n";
	}

private:

	bool error1(const PsimagLite::String& message, SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());
	}

	OptionsReadableType* optsReadable_;

}; // class InputCheck
} // namespace LanczosPlusPlus

/*@}*/
#endif
