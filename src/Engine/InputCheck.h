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
		} else if (label=="Connectors") {
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
			\item[TimeStepTargetting] TDMRG algorithm
			\item[DynamicTargetting] TBW
			\item[AdaptiveDynamicTargetting] TBW
			\item[CorrectionVectorTargetting] TBW
			\item[CorrectionTargetting] TBW
			\item[TargetingAncilla] TBW
			\item[MettsTargetting] TBW
			\item[TargetingInSitu] TBW
			\item[geometryallinsystem] During infinite algorithm make environment
			contain always exactly one site
			\item[vectorwithoffsets] TBW
			\item[allPvectors] TBW
			\item[printgeometry] TBW
			\item[tarEnable] Tar output files
			\item[tarNoDelete] Do not delete output files after taring them.
			\item[recoveryNoDelete] Do not delete recovery files even if run finishes OK
			\item[neverNormalizeVectors] TBW
			\item [advanceOnlyAtBorder] Advance time only at borders
			\item [findSymmetrySector] Find symmetry sector with lowest energy, and
			ignore value set in TargetElectronsUp or TargetSzPlusConst
			\item [KroneckerDumper] TBW
			\item [extendedPrint] TBW
			\item [useSvd] TBW
			\item [KronNoLoadBalance] Disable load balancing for MatrixVectorKron
			\item [setAffinities] TBW
			\item [wftInPatches] WFT calculation will be done using symmetry patches
			\item [diskstacks] Save and load stacks used for shrinking to and from disk,
							   instead of to and from memory
			\item [wftInBlocks] Accelerate the WFT by using dense blocks
			\item [wftStacksInDisk] Save and load stacks for WFT to and from disk,
							   instead of to and from memory. Cannot be used with restart yet.
			\item [BatchedGemm] Only meaningful with MatrixVectorKron. Enables
								batched gemm and might need plugin sc
			\item [KrylovAbridge] TBW
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
		registerOpts.push_back("TimeStepTargetting");
		registerOpts.push_back("DynamicTargetting");
		registerOpts.push_back("AdaptiveDynamicTargetting");
		registerOpts.push_back("CorrectionVectorTargetting");
		registerOpts.push_back("CorrectionTargetting");
		registerOpts.push_back("TargetingInSitu");
		registerOpts.push_back("TargetingRixsStatic");
		registerOpts.push_back("TargetingRixsDynamic");
		registerOpts.push_back("MettsTargetting");
		registerOpts.push_back("TargetingAncilla");
		registerOpts.push_back("geometryallinsystem");
		registerOpts.push_back("vectorwithoffsets");
		registerOpts.push_back("allPvectors");
		registerOpts.push_back("printgeometry");
		registerOpts.push_back("tarEnable");
		registerOpts.push_back("tarNoDelete");
		registerOpts.push_back("tarCoutNoDelete");
		registerOpts.push_back("recoveryNoDelete");
		registerOpts.push_back("normalizeTimeVectors");
		registerOpts.push_back("neverNormalizeVectors");
		registerOpts.push_back("noSaveStacks");
		registerOpts.push_back("noSaveData");
		registerOpts.push_back("noSaveWft");
		registerOpts.push_back("minimizeDisk");
		registerOpts.push_back("advanceOnlyAtBorder");
		registerOpts.push_back("findSymmetrySector");
		registerOpts.push_back("KroneckerDumper");
		registerOpts.push_back("doNotCheckTwoSiteDmrg");
		registerOpts.push_back("extendedPrint");
		registerOpts.push_back("useSvd");
		registerOpts.push_back("KronNoLoadBalance");
		registerOpts.push_back("setAffinities");
		registerOpts.push_back("wftInPatches");
		registerOpts.push_back("wftInBlocks");
		registerOpts.push_back("diskstacks");
		registerOpts.push_back("wftWithTemp");
		registerOpts.push_back("wftStacksInDisk");
		registerOpts.push_back("BatchedGemm");
		registerOpts.push_back("KrylovAbridge");

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

	PsimagLite::String getTargeting(const PsimagLite::String& options) const
	{
		PsimagLite::String targetting="GroundStateTargetting";

		const char *targets[]={"GroundStateTargetting",
		                       "TimeStepTargetting",
		                       "AdaptiveDynamicTargetting",
		                       "DynamicTargetting",
		                       "CorrectionVectorTargetting",
		                       "CorrectionTargetting",
		                       "MettsTargetting",
		                       "TargetingAncilla",
		                       "TargetingCorrelations",
		                       "TargetingInSitu",
		                       "TargetingRixsStatic",
		                       "TargetingRixsDynamic"};

		SizeType totalTargets = 12;

		SizeType count = 0;
		for (SizeType i = 0;i<totalTargets;++i) {
			if (options.find(targets[i])!=PsimagLite::String::npos) {
				if (targetting == "AdaptiveDynamicTargetting" &&
				        std::string(targets[i]) == "DynamicTargetting") continue;
				targetting = targets[i];
				count++;
			}
		}

		if (count > 1) {
			throw PsimagLite::RuntimeError("Only one targeting supported\n");
		}

		return targetting;
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

