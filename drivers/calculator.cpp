#include "ExpressionCalculator.h"

#ifdef USE_COMPLEX
typedef std::complex<double> ComplexOrRealType;
#else
typedef double ComplexOrRealType;
#endif

int main(int argc, char **argv)
{
	if (argc < 2) return 1;

	typedef PsimagLite::ExpressionCalculator<ComplexOrRealType> ExpressionCalculatorType;
	typedef PsimagLite::PrepassData<ComplexOrRealType> PrepassDataType;

	ExpressionCalculatorType::VectorStringType ve;
	PsimagLite::split(ve, argv[1], ",");


	PrepassDataType::VectorType vr(1,0.25);
	PrepassDataType pd("t", vr);

	PsimagLite::ExpressionPrepass<PrepassDataType>::prepass(ve,pd);

	ExpressionCalculatorType ec(ve);
	std::cout<<argv[1]<<"\t"<<ec()<<"\n";
}

