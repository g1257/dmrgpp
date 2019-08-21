/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
/** \ingroup DMRG */
/*@{*/

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "Options.h"
#include "Geometry/Geometry.h"
#include "ProgramGlobals.h"

namespace Dmrg {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	InputCheck() : optsReadable_(0)
	{
		allowedFileOptions_.push_back("");
		allowedFileOptions_.push_back("DELETE");
		allowedFileOptions_.push_back("list");
		allowedFileOptions_.push_back("keep");

		knownLabels_.push_back("TotalNumberOfSites");
		knownLabels_.push_back("NumberOfTerms");
		knownLabels_.push_back("DegreesOfFreedom");
		knownLabels_.push_back("GeometryKind");
		knownLabels_.push_back("GeometryOptions");
		knownLabels_.push_back("Connectors");
		knownLabels_.push_back("hubbardU");
		knownLabels_.push_back("potentialV");
		knownLabels_.push_back("Model");
		knownLabels_.push_back("SolverOptions");
		knownLabels_.push_back("Version");
		knownLabels_.push_back("OutputFile");
		knownLabels_.push_back("InfiniteLoopKeptStates");
		knownLabels_.push_back("FiniteLoops");
		knownLabels_.push_back("TargetElectronsUp");
		knownLabels_.push_back("TargetElectronsDown");
		knownLabels_.push_back("TargetSpinTimesTwo");
		knownLabels_.push_back("UseSu2Symmetry");
		knownLabels_.push_back("GsWeight");
		knownLabels_.push_back("TSPTau");
		knownLabels_.push_back("TSPTimeSteps");
		knownLabels_.push_back("TSPAdvanceEach");
		knownLabels_.push_back("ChebyshevTransform");
		knownLabels_.push_back("TSPAlgorithm");
		knownLabels_.push_back("TSPSites");
		knownLabels_.push_back("TSPLoops");
		knownLabels_.push_back("TSPProductOrSum");
		knownLabels_.push_back("IsPeriodicX");
		knownLabels_.push_back("Orbitals");
		knownLabels_.push_back("TSPOperator");
		knownLabels_.push_back("RAW_MATRIX");
		knownLabels_.push_back("FERMIONSIGN");
		knownLabels_.push_back("JMVALUES");
		knownLabels_.push_back("AngularFactor");
		knownLabels_.push_back("Threads");
		knownLabels_.push_back("LadderLeg");
		knownLabels_.push_back("linSize");
		knownLabels_.push_back("jvalues");
		knownLabels_.push_back("options");
		knownLabels_.push_back("version");
		knownLabels_.push_back("outputfile");
		knownLabels_.push_back("QNS");
		knownLabels_.push_back("BathSitesPerSite");
		knownLabels_.push_back("hoppings");
		knownLabels_.push_back("density");
		knownLabels_.push_back("skip");
		knownLabels_.push_back("ConnectorsX");
		knownLabels_.push_back("ConnectorsY");
		knownLabels_.push_back("ConnectorsBath");
		knownLabels_.push_back("RepeatFiniteLoopsTimes");
		knownLabels_.push_back("BetaDividedByTwo");
		knownLabels_.push_back("TSPRngSeed");
		knownLabels_.push_back("TSPOperatorMultiplier");
		knownLabels_.push_back("MettsCollapse");
		knownLabels_.push_back("HeisenbergTwiceS");
		knownLabels_.push_back("TargetElectronsTotal");
		knownLabels_.push_back("TargetSzPlusConst");
		knownLabels_.push_back("RepeatFiniteLoopsFrom");
		knownLabels_.push_back("RestartFilename");
		knownLabels_.push_back("RestartLabelForEnergy");
		knownLabels_.push_back("RestartMapStages");
		knownLabels_.push_back("RestartSourceTvForPsi");
		knownLabels_.push_back("RestartMappingTvs");
		knownLabels_.push_back("COOKED_OPERATOR");
		knownLabels_.push_back("COOKED_EXTRA");
		knownLabels_.push_back("OperatorExpression");
		knownLabels_.push_back("TargetExtra");
		knownLabels_.push_back("TSPEnergyForExp");
		knownLabels_.push_back("AdjustQuantumNumbers");
		knownLabels_.push_back("FeAsMode");
		knownLabels_.push_back("MaxMatrixRankStored");
		knownLabels_.push_back("CorrectionA");
		knownLabels_.push_back("DynamicDmrgType");
		knownLabels_.push_back("CorrectionVectorFreqType");
		knownLabels_.push_back("DynamicDmrgSteps");
		knownLabels_.push_back("DynamicDmrgEps");
		knownLabels_.push_back("DynamicDmrgAdvanceEach");
		knownLabels_.push_back("CorrectionVectorOmega");
		knownLabels_.push_back("CorrectionVectorEta");
		knownLabels_.push_back("CorrectionVectorAlgorithm");
		knownLabels_.push_back("CorrelationsType");
		knownLabels_.push_back("LongChainDistance");
		knownLabels_.push_back("IsPeriodicY");
		knownLabels_.push_back("TruncationTolerance");
		knownLabels_.push_back("PotentialT");
		knownLabels_.push_back("omega");
		knownLabels_.push_back("MagneticField");
		knownLabels_.push_back("SpinOrbit");
		knownLabels_.push_back("DegeneracyMax");
		knownLabels_.push_back("JzSymmetry");
		knownLabels_.push_back("DegeneracyMax");
		knownLabels_.push_back("KroneckerDumperBegin");
		knownLabels_.push_back("KroneckerDumperEnd");
		knownLabels_.push_back("LanczosEps");
		knownLabels_.push_back("LanczosSteps");
		knownLabels_.push_back("TridiagEps");
		knownLabels_.push_back("TridiagSteps");
		knownLabels_.push_back("TruncationTolerance");
		knownLabels_.push_back("GeometryMaxConnections");
		knownLabels_.push_back("LanczosNoSaveLanczosVectors");
		knownLabels_.push_back("DenseSparseThreshold");
		knownLabels_.push_back("TridiagonalEps");
		knownLabels_.push_back("HoneycombLy");
		knownLabels_.push_back("GeometryValueModifier");
		knownLabels_.push_back("ThreadsStackSize");
		knownLabels_.push_back("RecoverySave");
		knownLabels_.push_back("Intent");
		knownLabels_.push_back("PrintHamiltonianAverage");
		for (SizeType i = 0; i < 10; ++i)
			knownLabels_.push_back("Term" + ttos(i));
	}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		PsimagLite::String str = PsimagLite::Geometry<int,int,ProgramGlobals>::import();

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
		return str;
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		if (label=="JMVALUES" || label=="RS:JMVALUES") {
			if (vec.size()!=3) return error1("JMVALUES",line);
			return true;
		} else if (label=="RAW_MATRIX" || label=="RS:RAW_MATRIX" || label == "SpinOrbit") {
			if (!checkForMatrix(vec)) return error1(label,line);
			return true;
		} else if (label == "Connectors" || label == "hopOnSite") {
			if (!checkForMatrix(vec) && !checkForVector(vec))
				return error1(label,line);
			return true;
		} else if (label == "MagneticField") {
			return true;
		} else if (label=="FiniteLoops") {
			SizeType n = atoi(vec[0].c_str());
			if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
			return true;
		}

		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label,
	                      SizeType line) const
	{
		for (SizeType i = 0; i < knownLabels_.size(); ++i)
			if (knownLabels_[i] == label) return true;
		PsimagLite::String msg("WARNING: Unknown label " + label +"\n");
		std::cout<<msg;
		std::cerr<<msg;
		return false;
	}

	/* PSIDOC dmrgSolverOptions
	   \verb!SolverOptions=! in the input file must contain
		  a comma-separated list of strings. At least one of the following strings must
		  be provided:
		\begin{itemize}
			\item[none]  Use this when no options are given, because the list of
		   strings must be non-null.
				Note that ``none'' does not disable other options.

			 \item[useSu2Symmetry] Use the SU(2) symmetry for the model, and
			interpret quantum
				 numbers in the line ``QNS'' appropriately.

			 \item[nofiniteloops]  Don't do finite loops, even if provided under
			``FiniteLoops'' below.
			\item[restart] Restart from a previously saved run. See FIXME
			\item[debugmatrix] Print Hamiltonian matrix for targeted sector of
			superblock
			\item[exactdiag] Do exact diagonalization with LAPACK instead of Lanczos
			\item[nodmrgtransform] Do not DMRG transform bases
			\item[useDavidson] Use Davidson instead of Lanczos
			\item[verbose] Enable verbose output
			\item[nowft] Disable the Wave Function Transformation (WFT)
			\item[useComplex] TBW
			\item[inflate] TBW
			\item[twositedmrg] Use 2-site DMRG. Default is 1-site DMRG
			\item[noloadwft] TBW
			\item[ChebyshevSolver] Use ChebyshevSolver instead of Lanczos
			\item[MatrixVectorStored] Store superblock sector of Hamiltonian matrix
			in memory instead of constructing it on the fly.
			\item[MatrixVectorKron] TBW
			\item[TimeStepTargeting] TDMRG algorithm
			\item[DynamicTargeting] TBW
			\item[AdaptiveDynamicTargeting] TBW
			\item[CorrectionVectorTargeting] TBW
			\item[CorrectionTargeting] TBW
			\item[TargetingAncilla] TBW
			\item[MettsTargeting] TBW
			\item[TargetingInSitu] TBW
			\item[geometryallinsystem] During infinite algorithm make environment
			contain always exactly one site
			\item[vectorwithoffsets] TBW
			\item[allPvectors] TBW
			\item[printgeometry] TBW
			\item[recoveryEnableRead] Enables recovery if previous run crashed
			\item[neverNormalizeVectors] TBW
			\item [advanceUnrestricted] Don't restrict advance time to borders
			\item [findSymmetrySector] Find symmetry sector with lowest energy, and
			ignore value set in TargetElectronsUp or TargetSzPlusConst
			\item [KroneckerDumper] TBW
			\item [extendedPrint] TBW
			\item [truncationNoSvd] Do not use SVD for truncation;
									   use density matrix instead
			\item [KronNoLoadBalance] Disable load balancing for MatrixVectorKron
			\item [setAffinities] TBW
			\item [wftNoAccel] Disable WFT acceleration (but not the WFT itself)
			\item [wftAccelPatches] Force WFT acceleration with patches, even
			in twositedmrg
			\item [BatchedGemm] Only meaningful with MatrixVectorKron. Enables
								batched gemm and might need plugin sc
			\item [KrylovNoAbridge] TBW
			\item [fixLegacyBugs] TBW
			\item [saveDensityMatrixEigenvalues] Save DensityMatrixEigenvalues
												 to the data file.
			\item [KronNoUseLowerPart] Don't Use lower part of Kron matrix but
 recompute it instead.
			\item [shrinkStacksOnDisk] Store shrink stacks on disk instead of in memory
			\item [OperatorsChangeAll] Do not hollow out operators but keep track of
			them for all sites. This is will use more RAM, but might be needed
			to target expressions.
			\item [calcAndPrintEntropies] Calculate entropies and print to cout file
			\item [blasNotThreadSafe] TBW
		\end{itemize}
		*/
	void check(const PsimagLite::String& label,
	           const PsimagLite::String& val,
	           SizeType)
	{
		if (label!="SolverOptions") return;
		PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

		registerOpts.push_back("restart");
		registerOpts.push_back("debugmatrix");
		registerOpts.push_back("test");
		registerOpts.push_back("exactdiag");
		registerOpts.push_back("nodmrgtransform");
		registerOpts.push_back("useDavidson");
		registerOpts.push_back("verbose");
		registerOpts.push_back("nofiniteloops");
		registerOpts.push_back("nowft");
		registerOpts.push_back("useComplex");
		registerOpts.push_back("inflate");
		registerOpts.push_back("none");
		registerOpts.push_back("twositedmrg");
		registerOpts.push_back("noloadwft");
		registerOpts.push_back("ChebyshevSolver");
		registerOpts.push_back("MatrixVectorStored");
		registerOpts.push_back("MatrixVectorOnTheFly");
		registerOpts.push_back("TimeStepTargeting");
		registerOpts.push_back("DynamicTargeting");
		registerOpts.push_back("AdaptiveDynamicTargeting");
		registerOpts.push_back("CorrectionVectorTargeting");
		registerOpts.push_back("CorrectionTargeting");
		registerOpts.push_back("ChebyshevTargeting");
		registerOpts.push_back("TargetingInSitu");
		registerOpts.push_back("TargetingRixsStatic");
		registerOpts.push_back("TargetingRixsDynamic");
		registerOpts.push_back("MettsTargeting");
		registerOpts.push_back("TargetingAncilla");
		registerOpts.push_back("geometryallinsystem");
		registerOpts.push_back("vectorwithoffsets");
		registerOpts.push_back("allPvectors");
		registerOpts.push_back("printgeometry");
		registerOpts.push_back("recoveryEnableRead");
		registerOpts.push_back("normalizeTimeVectors");
		registerOpts.push_back("neverNormalizeVectors");
		registerOpts.push_back("noSaveStacks");
		registerOpts.push_back("noSaveData");
		registerOpts.push_back("noSaveWft");
		registerOpts.push_back("minimizeDisk");
		registerOpts.push_back("advanceUnrestricted");
		registerOpts.push_back("findSymmetrySector");
		registerOpts.push_back("KroneckerDumper");
		registerOpts.push_back("doNotCheckTwoSiteDmrg");
		registerOpts.push_back("extendedPrint");
		registerOpts.push_back("truncationNoSvd");
		registerOpts.push_back("KronNoLoadBalance");
		registerOpts.push_back("setAffinities");
		registerOpts.push_back("wftNoAccel");
		registerOpts.push_back("wftAccelPatches");
		registerOpts.push_back("BatchedGemm");
		registerOpts.push_back("KrylovNoAbridge");
		registerOpts.push_back("fixLegacyBugs");
		registerOpts.push_back("saveDensityMatrixEigenvalues");
		registerOpts.push_back("KronNoUseLowerPart");
		registerOpts.push_back("shrinkStacksOnDisk");
		registerOpts.push_back("OperatorsChangeAll");
		registerOpts.push_back("calcAndPrintEntropies");
		registerOpts.push_back("blasNotThreadSafe");

		PsimagLite::Options::Writeable optWriteable(registerOpts,
		                                            PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new  OptionsReadableType(optWriteable,val);

		bool mvs =  (val.find("MatrixVectorStored") != PsimagLite::String::npos);
		bool mvo =  (val.find("MatrixVectorOnTheFly") != PsimagLite::String::npos);
		bool notMvk = (mvs || mvo);
		if (val.find("BatchedGemm") != PsimagLite::String::npos) {
			if (notMvk)
				err("FATAL: BatchedGemm only with MatrixVectorKron\n");
#ifndef PLUGIN_SC
			err("BatchedGemm needs -DPLUGIN_SC in Config.make\n");
#endif
		}
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

	void checkFileOptions(PsimagLite::String fileOption)
	{
		if (passesFileOptions(fileOption)) return;
		throw PsimagLite::RuntimeError("Option " + fileOption + " not understood\n");
	}

private:

	bool passesFileOptions(PsimagLite::String fileOption)
	{
		for (SizeType i = 0; i < allowedFileOptions_.size(); ++i)
			if (std::find(allowedFileOptions_.begin(),
			              allowedFileOptions_.end(),
			              fileOption) == allowedFileOptions_.end())
				return false;
		return true;
	}

	bool checkForVector(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() == 0) return false;
		SizeType n = atoi(vec[0].c_str());
		return (vec.size() == n+1);
	}

	bool checkForMatrix(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() < 2) return false;
		SizeType row = atoi(vec[0].c_str());
		SizeType col = atoi(vec[1].c_str());
		SizeType n = row*col;
		return (vec.size() == n+2);
	}

	bool error1(const PsimagLite::String& message,SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());

	}

	OptionsReadableType* optsReadable_;
	VectorStringType allowedFileOptions_;
	VectorStringType knownLabels_;
}; // class InputCheck
} // namespace Dmrg

/*@}*/
#endif

