#ifndef RestartStruct_H
#define RestartStruct_H

#include "Vector.h"
#include <iostream>
#include "Io/IoSerializerStub.h"

namespace Dmrg {

struct RestartStruct {

	RestartStruct()
	: filename(""),into("GroundState"),labelForPsi("PSI"),labelForEnergy("Energy")
	{}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType = sizeof(RestartStruct),
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label,
	               PsimagLite::IoSerializer& ioSerializer) const
	{
		PsimagLite::String root = label;
		ioSerializer.createGroup(root);
		ioSerializer.write(root + "/filename", filename);
		ioSerializer.write(root + "/into", into);
		ioSerializer.write(root + "/labelForPsi", labelForPsi);
		ioSerializer.write(root + "/labelForEnergy", labelForEnergy);

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

