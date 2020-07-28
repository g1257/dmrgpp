/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef DMFT_INPUT_CHECK_H
#define DMFT_INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "../../PsimagLite/src/Options.h"
#include "Geometry/Geometry.h"
//#include "ProgramGlobals.h"

namespace Dmft {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	InputCheck() : optsReadable_(0)
	{
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
		knownLabels_.push_back("DynamicDmftType");
		knownLabels_.push_back("CorrectionVectorFreqType");
		knownLabels_.push_back("DynamicDmftSteps");
		knownLabels_.push_back("DynamicDmftEps");
		knownLabels_.push_back("DynamicDmftAdvanceEach");
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
		knownLabels_.push_back("SaveDensityMatrixEigenvalues");
	}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		//PsimagLite::String str = PsimagLite::Geometry<int,int,ProgramGlobals>::import();
		PsimagLite::String str;
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
		return false;
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::String& vec,
	           SizeType line) const
	{
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

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

private:

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
} // namespace Dmft

/*@}*/
#endif

