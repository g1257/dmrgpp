#ifndef _AINUR_EMPTY_H_
#define _AINUR_EMPTY_H_
#include <iostream>
#include <fstream>
#include "Vector.h"
#include "TypeToString.h"
#include "PsimagLite.h"
//#include "AinurStatements.h"

namespace PsimagLite {

class Ainur {

public:

	typedef Vector<String>::Type VectorStringType;
	//typedef AinurStatements AinurStatementsType;
	//typedef AinurStatementsType::AinurLexicalType AinurLexicalType;

	Ainur(String str)
	    : dummy_("")
	{
		errorMessage();
	}

	String& prefix() { return dummy_; }

	const String& prefix() const { return dummy_; }

	void printUnused(std::ostream& os) const
	{
		errorMessage();
	}

	template<typename SomeType>
	void readValue(SomeType& t, String label) const
	{
		errorMessage();
	}

private:

	void errorMessage() const
	{
		err("To use Ainur, you need boost-devel, and compile with -DUSE_BOOST\n");
	}

	String dummy_;
};

}
#endif // _AINUR_EMPTY_H_
