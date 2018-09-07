#include "LanczosSolver.h"
#include "Matrix.h"
#include "CrsMatrix.h"
#include <fstream>

int main(int argc, char **argv)
{
	if (argc != 2) {
		std::cerr<<"Expected filename\n";
		return 1;
	}

	typedef double RealType;
	typedef std::complex<RealType> ComplexType;
	typedef PsimagLite::ParametersForSolver<double> SolverParametersType;
	typedef PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorType;

	PsimagLite::Matrix<ComplexType> m;
	std::ifstream fin(argv[1]);
	PsimagLite::String file(argv[1]);
	if (!fin || fin.bad() || !fin.good())
		throw PsimagLite::RuntimeError("Could not open file " + file + "\n");

	SizeType tmp = 0;
	fin>>tmp;
	SizeType n = 1;
	fin>>n;
	if (tmp != n)
		throw PsimagLite::RuntimeError("Matrix in file " + file + " is not square\n");

	m.resize(n, n , 0);
	for (SizeType i = 0; i < n; ++i)
		for (SizeType j = 0; j < n; ++j)
			fin>>m(i, j);

	PsimagLite::CrsMatrix<ComplexType> msparse(m);

	SolverParametersType params;
	params.lotaMemory = true;
	params.tolerance = -1;
	params.options = "reortho";

	PsimagLite::LanczosSolver<SolverParametersType,
	        PsimagLite::CrsMatrix<ComplexType>,VectorType> lanczosSolver(msparse, params);

	VectorType z(n, 0.0);
	VectorRealType e(n);
	for (SizeType excited = 0; excited < n; ++excited)
		lanczosSolver.computeExcitedState(e[excited], z, excited);

	VectorRealType eigs(n);
	PsimagLite::diag(m, eigs, 'V');
	std::cout<<"\nEXACT: ";
	for (SizeType excited = 0; excited < n; ++excited)
		std::cout<<eigs[excited]<<" ";
	std::cout<<"\n";

	std::cout<<"LANCZOS: ";
	for (SizeType excited = 0; excited < n; ++excited)
		std::cout<<e[excited]<<" ";
	std::cout<<"\n";
}
