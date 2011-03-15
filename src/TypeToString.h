
/* PsimagLite */
/* See LICENSE for licensing details and disclaimers */
#ifndef TYPE_TO_STRING_H
#define TYPE_TO_STRING_H

#include <string>
#include <sstream>

namespace PsimagLite {
	template<class T>
	std::string typeToString(T t)
	{
		std::stringstream ss;
		std::string str;
		ss.precision(10);
		ss<<t;
		ss>>str;
		return str;
	}
}

template<class T>
std::string ttos(T t)
{
	return PsimagLite::typeToString(t);
}

#endif // TYPE_TO_STRING_H

