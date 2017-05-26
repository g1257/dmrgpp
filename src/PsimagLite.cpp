#include "PsimagLite.h"

namespace PsimagLite {

std::ostream& operator<<(std::ostream& os,const std::pair<SizeType,SizeType>& p)
{
	os<<p.first<<" "<<p.second<<" ";
	return os;
}

std::istream& operator>>(std::istream& is,std::pair<SizeType,SizeType>& pair)
{
	is>>pair.first;
	is>>pair.second;
	return is;
}

void err(String s)
{
	throw RuntimeError(s);
}

SizeType log2Integer(SizeType x)
{
	SizeType count = 0;
	while (x > 0) {
		x >>= 1;
		++count;
	}

	return count;
}

String basename(const String& path)
{
	return String(std::find_if(path.rbegin(),
	                           path.rend(),
	                           MatchPathSeparator()).base(),
	              path.end());
}

const int PsiApp::libSizeOfSizeType_ = sizeof(SizeType);

} // namespace PsimagLite

void err(PsimagLite::String s)
{
	PsimagLite::err(s);
}

