#include "LanczosSolver.h"
#include "Matrix.h"
#include "CrsMatrix.h"
#include <fstream>
#include "PsimagLite.h"

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

	typedef PsimagLite::LanczosSolver<SolverParametersType,
	        PsimagLite::CrsMatrix<ComplexType>,VectorType> LanczosSolverType;

	LanczosSolverType lanczosSolver(msparse, params);

	VectorType initial(n);
	PsimagLite::fillRandom(initial);

	VectorRealType eigs(n);
	PsimagLite::diag(m, eigs, 'V');
	std::cout<<"\nEXACT: ";
	for (SizeType excited = 0; excited < n; ++excited)
		std::cout<<eigs[excited]<<" ";
	std::cout<<"\n";

	std::cout<<"LANCZ: ";
	LanczosSolverType::VectorVectorType zz;
	lanczosSolver.computeAllStatesBelow(eigs, zz, initial, n);

	std::cout<<"LANCZOS: \n";
	for (SizeType excited = 0; excited < n; ++excited) {
		std::cout<< "<E>_" << excited << " = " << eigs[excited] <<"  \n";
	}

	std::cout<<"\n";
}
