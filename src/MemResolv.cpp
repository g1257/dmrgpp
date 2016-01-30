#include "MemResolv.h"

namespace PsimagLite {

bool operator==(const MemoryPointer& a, const MemoryPointer& b)
{
	bool b1 = (a.type == b.type);
	bool b2 = (a.length == b.length);
	bool b3 = (a.ptr == b.ptr);

	return (b1 & b2 & b3);
}

std::ostream& operator<<(std::ostream& os, const MemResolv& mresolv)
{
	os<<"MemResolvs "<<mresolv.vmptr_.size()<<"\n";
	for (SizeType i = 0; i < mresolv.vmptr_.size(); ++i)
		mresolv.print(os,mresolv.vmptr_[i]);

	os<<"MemResolv garbage: "<<mresolv.garbage_.size();
	for (SizeType i = 0; i < mresolv.garbage_.size(); ++i) {
		os<<reinterpret_cast<void *>(mresolv.garbage_[i]);
		os<<" "<<mresolv.garbageSize_[i];
	}

	os<<"\n";

	return os;
}

} // namespace PsimagLite

