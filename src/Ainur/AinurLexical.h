#ifndef AINURLEXICAL_H
#define AINURLEXICAL_H
#include "Vector.h"

namespace PsimagLite {

class AinurLexical {

public:

	typedef Vector<String>::Type VectorStringType;

	static void removeTrailingWhitespace(VectorStringType& v)
	{
		for (SizeType i = 0; i < v.size(); ++i)
			removeTrailingBlanks(v[i]);
	}

	static void removeTrailingBlanks(String& s)
	{
		SizeType start = 0;
		SizeType l = s.length();
		for (SizeType i = 0; i < l; ++i) {
			if (isWhitespace(s[i]) || isEOL(s[i]))
				start = i + 1;
			else
				break;
		}

		if (start == l) {
			s = "";
			return;
		}

		String newStr = s.substr(start);
		l = newStr.length();
		SizeType end = l;
		for (int i = l - 1; i >= 0; --i) {
			if (isWhitespace(newStr[i]) || isEOL(newStr[i]))
				end = i;
			else
				break;
		}

		s = newStr.substr(0, end);
	}

	static bool isWhitespace(char c)
	{
		return (c == ' ' || c == '\t');
	}

	static bool isEOL(char c)
	{
		return (c == '\n' || c == '\r');
	}
}; // class AinurLexical
} // namespace PsimagLite
#endif // AINURLEXICAL_H
