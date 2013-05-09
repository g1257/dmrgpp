
/* PsimagLite */
/* See LICENSE for licensing details and disclaimers */
#ifndef TYPE_TO_STRING_H
#define TYPE_TO_STRING_H

#include <sstream>
#include "String.h"

namespace PsimagLite {
	template<class T>
	String typeToString(T t)
	{
		std::stringstream ss;
		String str;
		ss.precision(10);
		ss<<t;
		ss>>str;
		return str;
	}
}

template<class T>
PsimagLite::String ttos(T t)
{
	return PsimagLite::typeToString(t);
}

#endif // TYPE_TO_STRING_H

