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

void split(Vector<String>::Type& tokens, String str, String delimiters)
{
	// Skip delimiters at beginning.
	String::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	String::size_type pos     = str.find_first_of(delimiters, lastPos);

	while (String::npos != pos || String::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

String basename(const String& path)
{
	return String(std::find_if(path.rbegin(),
	                           path.rend(),
	                           MatchPathSeparator()).base(),
	              path.end());
}

int atoi(String str)
{
	const SizeType n = str.length();
	for (SizeType i = 0; i < n; ++i) {
		if (isdigit(str[i])) continue;
		if (str[i] == '-' || str[i] == '+') continue;
		throw RuntimeError("atoi received a non-digit\n");
	}

	return std::atoi(str.c_str());
}

bool isAnInteger(String str)
{
	const SizeType n = str.length();
	for (SizeType i = 0; i < n; ++i) {
		if (isdigit(str[i])) continue;
		return false;
	}

	return true;
}

bool isAfloat(String str)
{
	const SizeType n = str.length();
	for (SizeType i = 0; i < n; ++i) {
		if (isdigit(str[i])) continue;
		if (str[i] == '-' ||
		        str[i] == '+' ||
		        str[i] == 'e' ||
		        str[i] == 'E' ||
		        str[i] == '.') continue;
		return false;
	}

	return true;
}

double atof(String str)
{

	if (!isAfloat(str))
		throw RuntimeError("atof received a non-digit " + str + "\n");

	return std::atof(str.c_str());
}

const int PsiApp::libSizeOfSizeType_ = sizeof(SizeType);

} // namespace PsimagLite

void err(PsimagLite::String s)
{
	PsimagLite::err(s);
}

