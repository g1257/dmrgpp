#include <cassert>
#include <iostream>
#include <sstream>

using RealType = double;
using SizeType = unsigned int;

RealType calcMu(SizeType site, SizeType n, SizeType N1, RealType tau, RealType mu)
{
	//  (nm is the number of steps to increase the onsite
	//   chemical potential at site i)
	static SizeType nm = 0;
	static SizeType rm = 1; //( steps changes when n=1000)

	if (site + 1 == N1 || site == 0) {
		++nm;

		if (n % 1000 == 0 && n > 0) {
			++rm;
			nm = 0;

			std::stringstream msg;
			msg << "FermionSpinless::calcMu(): nm=0 "
			    << "site=" << site << " n=" << n << " tau=" << tau << " mu=" << mu
			    << " rm=" << rm << "\n";
			std::cout << msg.str();
		}
	}

	if (site + rm < N1) { //(chemical potential for  i <= end site - rm)
		return -mu;
	}

	if (site + rm > N1) { //(chemical potential for  i > end site -rm+1)
		return -64.0;
	}

	assert(site + rm == N1); //(chemical potential for  i = end site -rm+1)
	return -32.0 * tau * nm;
}

int main(int argc, char* argv[])
{
	RealType tau = 0.001;
	RealType mu = 0;
	SizeType nsites = 12;
	SizeType totalTimeSteps = 2000;

	for (SizeType currentTimeStep = 0; currentTimeStep < totalTimeSteps; ++currentTimeStep) {
		for (SizeType site = 0; site < nsites; ++site) {
			const RealType effectiveMu = calcMu(site, currentTimeStep, nsites, tau, mu);
			std::cout << effectiveMu << " ";
		}

		RealType time = tau * currentTimeStep;
		std::cout << ": TIME=" << time << "\n";
	}
}
