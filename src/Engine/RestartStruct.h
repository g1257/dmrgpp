#ifndef RestartStruct_H
#define RestartStruct_H

#include "Vector.h"
#include <iostream>

namespace Dmrg {

struct RestartStruct {

	RestartStruct() : filename(""),into("GroundState"),labelForPsi("PSI")
	{}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType = sizeof(RestartStruct),
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	PsimagLite::String filename;
	PsimagLite::String into;
	PsimagLite::String labelForPsi;
};

std::ostream& operator<<(std::ostream& os, const RestartStruct& c);

std::istream& operator>>(std::istream& is,RestartStruct& c);

} // namespace Dmrg

#endif // RestartStruct_H

