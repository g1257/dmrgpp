#ifndef RestartStruct_H
#define RestartStruct_H

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

} // namespace Dmrg

#endif // RestartStruct_H

