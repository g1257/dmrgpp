#include "AST/PlusMinusMultiplyDivide.h"
#include "AST/ExpressionForAST.h"
#include "PredicateAwesome.h"
#include "PsimagLite.h"

#ifdef USE_COMPLEX
typedef std::complex<double> ComplexOrRealType;
#else
typedef double ComplexOrRealType;
#endif

int main(int argc, char **argv)
{
	if (argc < 2) return 1;

	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::PlusMinusMultiplyDivide<ComplexOrRealType> PrimitivesType;
	PsimagLite::String str(argv[1]);
	PsimagLite::replaceAll(str, "%t", "0.25");
	VectorStringType ve;
	PsimagLite::split(ve, str, ":");

	PrimitivesType primitives;
	PsimagLite::ExpressionForAST<PrimitivesType> expresionForAST(ve, primitives);

	std::cout<<argv[1]<<"\t"<<expresionForAST.exec()<<"\n";
}

