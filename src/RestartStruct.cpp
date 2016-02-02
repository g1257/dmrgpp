#include "RestartStruct.h"

namespace Dmrg {

std::ostream& operator<<(std::ostream& os, const RestartStruct& c)
{
	if (c.filename == "") return os;

	os<<"RestartStruct.filename="<<c.filename<<"\n";
	os<<"RestartStruct.into="<<c.into<<"\n";
	os<<"RestartStruct.labelForPsi="<<c.labelForPsi<<"\n";
	os<<"RestartStruct.labelForEnergy="<<c.labelForEnergy<<"\n";
	return os;
}

std::istream& operator>>(std::istream& is,RestartStruct& c)
{
	is>>c.filename;
	is>>c.into;
	is>>c.labelForPsi;
	is>>c.labelForEnergy;
	return is;
}

} // namespace Dmrg

