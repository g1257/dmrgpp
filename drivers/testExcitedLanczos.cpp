#include "LanczosSolver.h"
#include "Matrix.h"
#include "CrsMatrix.h"
#include "Random48.h" // <--- is this thread safe?

int main(int argc,char *argv[])
{
	typedef double RealType;
	typedef std::complex<double> ComplexType;
	typedef PsimagLite::ParametersForSolver<double> SolverParametersType;
	typedef PsimagLite::Vector<ComplexType>::Type VectorComplexType;
	typedef PsimagLite::Vector<RealType>::Type VectorRealType;

	if (argc!=2) {
		std::cout<<"USAGE: "<<argv[0]<<" excited\n";
		return 1;
	}

	SizeType excited = atoi(argv[1]);

	SizeType n = 4;
	PsimagLite::Matrix<ComplexType> m(n, n);
	// fill m
	m(1,2) = -0.5;
	m(2,1) = PsimagLite::conj(m(1,2));
	m(0,0) = m(3,3) = -0.25;
	m(1,1) = m(2,2) = 0.25;

	std::cout<<"matrix\n"<<m<<"\n";

	PsimagLite::CrsMatrix<ComplexType> msparse(m);

	SolverParametersType params;
	params.lotaMemory = true;
	params.tolerance = -1;
	params.options = "reortho";
	PsimagLite::LanczosSolver<SolverParametersType,
	        PsimagLite::CrsMatrix<ComplexType>,VectorComplexType> lanczosSolver(msparse, params);


	PsimagLite::Random48<double> myrng(time(0));
	VectorComplexType initialV(n, 0.0);
	for (SizeType i = 0; i < n; ++i)
		initialV[i] = myrng() - 0.5;

	double e1 = 0;
	VectorComplexType z1(n, 0.0);
	lanczosSolver.computeOneState(e1, z1, initialV, excited);
	std::cout<<"energy1="<<e1<<"\n";
	std::cout<<z1<<"\n";

	VectorRealType eigs(m.rows());
	PsimagLite::diag(m,eigs,'V');
	for (SizeType i=0;i<m.rows();i++) {
		std::cerr<<"eigs["<<i<<"]="<<eigs[i]<<"\n";
	}
	VectorRealType zeig(n, 0.0);
	for (SizeType i=0;i<m.rows();i++) {
		zeig[i] =  PsimagLite::real(m(i,0));
	}
	std::cout<<m<<"\n";
}
