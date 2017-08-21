#include "Ainur.h"

PsimagLite::String dmrgImport()
{
	PsimagLite::String str("");
	str += "integer TotalNumberOfSites;\n";
	str += "integer NumberOfTerms;\n";
	str += "integer DegreesOfFreedom;\n";
	str += "string GeometryKind;\n";
	str += "string GeometryOptions;\n";
	str += "vector dir0:Connectors;\n";
	str += "vector dir1:Connectors;\n";
	str += "integer LadderLeg;\n";
	str += "vector hubbardU;\n";
	str += "vector potentialV;\n";
	str += "string Model;\n";
	str += "string SolverOptions;\n";
	str += "string Version;\n";
	str += "integer InfiniteLoopKeptStates;\n";
	str += "string OutputFile;\n";
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
	// str += "integer Threads = 1;\n"; //<--- NOT ALLOWED YET

	return str;
}

void partiallyReadSomething(const PsimagLite::Ainur& ainur)
{
	SizeType n = 0;
	ainur.readValue(n, "TotalNumberOfSites");
	std::cout<<"TotalNumberOfSites="<<n<<"\n";
}

int main(int argc, char** argv)
{
//	if (argc == 1) return 1;
//	std::ifstream fin(argv[1]);
	PsimagLite::String str;

//	fin.seekg(0, std::ios::end);
//	str.reserve(fin.tellg());
//	fin.seekg(0, std::ios::beg);

//	str.assign((std::istreambuf_iterator<char>(fin)),
//	           std::istreambuf_iterator<char>());
//	fin.close();

	str = dmrgImport() + str;
	PsimagLite::Ainur ainur(str);
	partiallyReadSomething(ainur);
	ainur.printAll(std::cout);
}
