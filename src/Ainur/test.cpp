#include "Ainur.h"

PsimagLite::String dmrgImport()
{
	PsimagLite::String str("");
	str += "require TotalNumberOfSites.typeof = \"integer\";\n";
	str += "require NumberOfTerms.typeof = \"integer\";\n";
	str += "require DegreesOfFreedom.typeof = \"integer\";\n";
	str += "require GeometryKind.typeof = \"string\";\n";
	str += "require GeometryOptions.typeof = \"string\";\n";
	str += "let dir0.Connectors;\n";
	str += "let dir1.Connectors;\n";
	str += "require LadderLeg.typeof = \"integer\";\n";
	str += "let hubbardU.typeof = \"vector\";\n";
	str += "let potentialV.typeof = \"vector\";\n";
	str += "require Model.typeof = \"string\";\n";
	str += "require SolverOptions.typeof = \"string\";\n";
	str += "require Version.typeof = \"string\";\n";
	str += "require InfiniteLoopKeptStates.typeof = \"integer\";\n";
	str += "require OutputFile.typeof = \"string\";\n";
	str += "require FiniteLoops.typeof = \"matrix\";\n";
	str += "require RepeatFiniteLoopsFrom.typeof = \"integer\";\n";
	str += "require RepeatFiniteLoopsTimes.typeof = \"integer\";\n";
	str += "let TargetElectronsUp.typeof = \"integer\";\n";
	str += "let TargetElectronsDown.typeof = \"integer\";\n";
	str += "let GsWeight.typeof = \"real\";\n";
	str += "let TSPTau.typeof = \"real\";\n";
	str += "let TSPTimeSteps.typeof = \"integer\";\n";
	str += "let TSPAdvanceEach.typeof = \"integer\";\n";
	str += "let TSPAlgorithm.typeof = \"string\";\n";
	str += "let TSPSites.typeof = \"vector.integer\";\n";
	str += "let TSPLoops.typeof = \"vector.integer\";\n";
	str += "let TSPProductOrSum.typeof = \"string\";\n";
	str += "let TSPOperator.typeof = \"string\";\n";
	str += "let OperatorExpression.typeof = \"string\";\n";
	str += "let Threads.typeof = \"integer\";\n";
	return str;
}

int main(int argc, char** argv)
{
	if (argc == 1) return 1;
	PsimagLite::String dmrgppImport = dmrgImport();
	PsimagLite::Ainur ainur(argv[1], dmrgppImport);
}
