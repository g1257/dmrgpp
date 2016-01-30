#include "FiniteLoop.h"

namespace Dmrg {

std::istream &operator>>(std::istream& is,FiniteLoop& fl)
{
	is>>fl.stepLength;
	is>>fl.keptStates;
	is>>fl.saveOption;
	return is;
}

std::ostream &operator<<(std::ostream& os,const FiniteLoop& fl)
{
	os<<fl.stepLength<<" ";
	os<<fl.keptStates<<" ";
	os<<fl.saveOption;
	return os;
}

} // namespace Dmrg

