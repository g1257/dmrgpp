#ifndef RestartStruct_H
#define RestartStruct_H

#include "Vector.h"
#include <iostream>
#include "IoSerializerStub.h"

namespace Dmrg {

struct RestartStruct {

	RestartStruct()
	: filename(""),into("GroundState"),labelForPsi("PSI"),labelForEnergy("#Energy")
	{}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType = sizeof(RestartStruct),
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	void serialize(PsimagLite::String label,
	               PsimagLite::IoSerializer& ioSerializer) const
	{
		PsimagLite::String root = label;
		ioSerializer.createGroup(root);
		ioSerializer.writeToTag(root + "/filename", filename);
		ioSerializer.writeToTag(root + "/into", into);
		ioSerializer.writeToTag(root + "/labelForPsi", labelForPsi);
		ioSerializer.writeToTag(root + "/labelForEnergy", labelForEnergy);

	}

	PsimagLite::String filename;
	PsimagLite::String into;
	PsimagLite::String labelForPsi;
	PsimagLite::String labelForEnergy;
};

std::ostream& operator<<(std::ostream& os, const RestartStruct& c);

std::istream& operator>>(std::istream& is,RestartStruct& c);

} // namespace Dmrg

#endif // RestartStruct_H

