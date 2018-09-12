#include "LanczosSolver.h"
#include "Matrix.h"
#include "CrsMatrix.h"
#include "PsimagLite.h"

int main()
{
	typedef std::complex<double> ComplexType;
	typedef PsimagLite::ParametersForSolver<double> SolverParametersType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorType;

	int n = 8;
	PsimagLite::Matrix<ComplexType> m(n, n);
	// fill m
	for (int i = 0; i < n; ++i) m(i,i) = 1.0;
	m(1,2) = ComplexType(0.0,0.5);
	m(2,1) = PsimagLite::conj(m(1,2));
	m(3,5) = m(5,3) = -1.5;

	PsimagLite::CrsMatrix<ComplexType> msparse(m);

	SolverParametersType params;
	params.lotaMemory = true;

	PsimagLite::LanczosSolver<SolverParametersType,
	        PsimagLite::CrsMatrix<ComplexType>,VectorType> lanczosSolver(msparse, params);

	double e = 0;
	VectorType z(n, 0.0);
	VectorType initial(n);
	PsimagLite::fillRandom(initial);
	lanczosSolver.computeOneState(e, z, initial, 0);

	std::cout<<"energy="<<e<<"\n";
	std::cout<<z<<"\n";
}
