#include "Ainur.h"

PsimagLite::String dmrgImport()
{
	PsimagLite::String str("");
	str += "require integer TotalNumberOfSites;\n";
	str += "require integer NumberOfTerms;\n";
	str += "require integer DegreesOfFreedom;\n";
	str += "require string GeometryKind;\n";
	str += "require string GeometryOptions;\n";
	str += "vector dir0:Connectors;\n";
	str += "vector dir1:Connectors;\n";
	str += "integer LadderLeg;\n";
	str += "vector hubbardU;\n";
	str += "vector potentialV;\n";
	str += "require string Model;\n";
	str += "require string SolverOptions;\n";
	str += "require string Version;\n";
	str += "require integer InfiniteLoopKeptStates;\n";
	str += "require string OutputFile;\n";
	str += "require matrix.integer FiniteLoops;\n";
	str += "require integer RepeatFiniteLoopsFrom;\n";
	str += "require integer RepeatFiniteLoopsTimes;\n";
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

	return str;
}

void partiallyReadSomething(const PsimagLite::Ainur& ainur)
{
	SizeType n = ainur.getLabel("TotalNumberOfSites");
	std::cout<<"TotalNumberOfSites="<<n<<"\n";
}

int main(int argc, char** argv)
{
	if (argc == 1) return 1;
	PsimagLite::String dmrgppImport = dmrgImport();
	PsimagLite::Ainur ainur(argv[1], dmrgppImport);
	partiallyReadSomething(ainur);

}
