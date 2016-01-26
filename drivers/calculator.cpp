#include "ExpressionCalculator.h"
#include "Tokenizer.h"

int main(int argc, char **argv)
{
	if (argc < 2) return 1;

	typedef PsimagLite::ExpressionCalculator<double> ExpressionCalculatorType;
	typedef PsimagLite::PrepassData<double> PrepassDataType;

	ExpressionCalculatorType::VectorStringType ve;
	PsimagLite::tokenizer(argv[1],ve,",");

	PrepassDataType pd;
	PrepassDataType::VectorType vr(1,0.25);
	pd.names = "t";
	pd.values = vr;

	PsimagLite::ExpressionPrepass<PrepassDataType>::prepass(ve,pd);

	ExpressionCalculatorType ec(ve);
	std::cout<<argv[1]<<"\t"<<ec()<<"\n";
}

