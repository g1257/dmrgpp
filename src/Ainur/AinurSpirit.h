#ifndef _AINUR_SPIRIT_H_
#define _AINUR_SPIRIT_H_
#include <iostream>
#include <fstream>
#include "Vector.h"
#include "TypeToString.h"
#include "PsimagLite.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <iostream>
#include <string>
#include <fstream>

namespace PsimagLite {

class Ainur {

public:

	typedef Vector<String>::Type VectorStringType;
	//typedef AinurStatements AinurStatementsType;
	//typedef AinurStatementsType::AinurLexicalType AinurLexicalType;

	Ainur(String str)
	    : dummy_("")
	{}

	String& prefix() { return dummy_; }

	const String& prefix() const { return dummy_; }

	void printUnused(std::ostream& os) const
	{
		os<<"PRINT UNUSED\n";
	}

	template<typename SomeType>
	void readValue(SomeType& t, String label) const
	{
		std::cerr<<"readValue called for label="<<label<<"\n";
	}

private:

	String dummy_;
};

}
#endif // _AINUR_SPIRIT_H_
