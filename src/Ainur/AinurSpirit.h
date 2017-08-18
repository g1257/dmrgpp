#ifndef _AINUR_SPIRIT_H_
#define _AINUR_SPIRIT_H_
#include <iostream>
#include <fstream>
#include "../Vector.h"
#include "../TypeToString.h"
#include "../PsimagLite.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
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

	typedef boost::fusion::vector<std::string, std::string> AttribType;
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
	void readValue(SomeType& t, String label) const;

private:

	String dummy_;
	boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::qi::space_type> keywords_;
	boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::qi::space_type> quoted2_;
	boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::qi::space_type> quotedString_;
	boost::spirit::qi::rule<Iterator, AttribType, boost::spirit::qi::space_type> statement1_;
}; // class AinurSpirit

}
#endif // _AINUR_SPIRIT_H_
