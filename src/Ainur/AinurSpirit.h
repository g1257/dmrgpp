#ifndef _AINUR_SPIRIT_H_
#define _AINUR_SPIRIT_H_
#include <iostream>
#include <fstream>
#include "../Vector.h"
#include "../TypeToString.h"
#include "../PsimagLite.h"
#include <iostream>
#include <string>
#include <fstream>

namespace PsimagLite {

class Ainur {

	struct Action {

		Action(const char *name)
		    : name_(name)
		{}

		template <typename A, typename ContextType, typename PassType>
		void operator()(A& attr,
		                ContextType& context,
		                PassType hit) const;

	private:

		const char* name_;
	};

public:

	typedef std::string::iterator Iterator;
	typedef Vector<String>::Type VectorStringType;
	//typedef AinurStatements AinurStatementsType;
	//typedef AinurStatementsType::AinurLexicalType AinurLexicalType;

	Ainur(String str);

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
}; // class AinurSpirit

}
#endif // _AINUR_SPIRIT_H_
