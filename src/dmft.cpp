#include "DmftSolver.h"
#include "Dispersion.h"

int main()
{
	typedef Dmft::DmftSolver<std::complex<double> > DmftSolverType;
	typedef DmftSolverType::RealType RealType;
	typedef DmftSolverType::DispersionType DispersionType;

	RealType ficticiousBeta = 0.01;
	SizeType nMatsubara = 100;
	DispersionType dispersion;

	DmftSolverType dmftSolver(ficticiousBeta, nMatsubara, dispersion);

	dmftSolver.selfConsistencyLoop();
}
