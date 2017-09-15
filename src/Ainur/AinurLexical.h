#ifndef AINURLEXICAL_H
#define AINURLEXICAL_H
#include "../Vector.h"

/* This class only checks whether input file contains valid characters
 * and throws if not */
namespace PsimagLite {

class AinurLexical {

public:

	AinurLexical(const String& str)
	{
		for (String::const_iterator it = str.begin(); it != str.end(); ++it) {
			if (allowedChar(*it)) continue;
			unsigned int x = *it;
			err("Lexical error: Not allowed char with code " +
			    ttos(x) + " near " + getContext(it,str.begin(),str.end()) + "\n");
		}
	}

	bool static isEmptyChar(char c)
	{
		return 	(isWhitespace(c) || isEOL(c));
	}

private:

	String getContext(String::const_iterator c,
	                  String::const_iterator b,
	                  String::const_iterator e,
	                  SizeType n = 10) const
	{
		String::const_iterator begin = (c - n > b) ? c -n : b;
		String::const_iterator end = (c + n > e) ? e : c + n;

		return String(begin, end);
	}

	static bool allowedChar(unsigned char c)
	{
		if (isWhitespace(c) || isEOL(c)) return true;
		return (c < 33 || c > 126 || c == 96) ? false : true;
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
