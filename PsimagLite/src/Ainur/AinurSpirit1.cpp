#include "AinurLexical.h"
#include "AinurSpirit.h"
#include <boost/config/warning_disable.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/mpl.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/spirit/include/qi.hpp>

namespace PsimagLite {

template <typename A, typename ContextType>
void Ainur::Action::operator()(A& attr, ContextType&, bool&) const
{
	if (name_ == "statement1") {
		String v1 = boost::fusion::at_c<0>(attr);
		String v2 = boost::fusion::at_c<1>(attr);
		state_.assign(v1, v2);
	} else if (name_ == "statement2") {
		String v1 = boost::fusion::at_c<0>(attr);
		String v2 = boost::fusion::at_c<1>(attr);
		state_.declare(v1, v2);
	} else {
		err("Ainur: bad action name " + name_ + "\n");
	}
}

template <typename A, typename ContextType>
void Ainur::Action3::operator()(A& attr, ContextType&, bool&) const
{
	String v1 = boost::fusion::at_c<0>(attr);
	String v2 = boost::fusion::at_c<1>(attr);
	String v3 = boost::fusion::at_c<2>(attr);
	state_.declare(v1, v2, v3);
}

Ainur::Ainur(String str)
    : dummy_("")
{
	AinurLexical discard(str);
	bool         verbose = AinurState::verbose();

	if (verbose)
		std::cerr << str << "\n\n";
#define AINUR_COMMENTS ('#' >> *(qi::char_ - qi::eol) >> qi::eol) | qi::eol | qi::space
	namespace qi      = boost::spirit::qi;
	namespace ascii   = boost::spirit::ascii;
	using AttribType  = boost::fusion::vector<std::string, std::string>;
	using Attrib3Type = boost::fusion::vector<std::string, std::string, std::string>;

	using SkipperType = BOOST_TYPEOF(AINUR_COMMENTS);

	qi::rule<IteratorType, std::string(), qi::unused_type> aToZ;
	qi::rule<IteratorType, std::string(), qi::unused_type> zeroToNine;
	qi::rule<IteratorType, std::string(), qi::unused_type> keywords;
	qi::rule<IteratorType, std::string(), qi::unused_type> value;
	qi::rule<IteratorType, std::string(), qi::unused_type> typeQualifier;
	qi::rule<IteratorType, AttribType, SkipperType>        statement1;
	qi::rule<IteratorType, AttribType, SkipperType>        statement2;
	qi::rule<IteratorType, Attrib3Type, SkipperType>       statement3;
	qi::rule<IteratorType, SkipperType>                    statement;
	value %= qi::lexeme[+(qi::char_ - qi::char_(";"))];
	aToZ       = ascii::char_("a", "z") | ascii::char_("A", "Z");
	zeroToNine = ascii::char_("0", "9");

	typeQualifier %= +(aToZ | ascii::char_(".") | ascii::char_("!"));
	keywords = *(ascii::char_("_")) >> +aToZ
	    >> *(ascii::char_("a", "z") | ascii::char_("A", "Z") | ascii::char_("0", "9")
	         | ascii::char_(":") | ascii::char_("_"));
	statement1 %= keywords >> '=' >> value;
	statement2 %= typeQualifier >> keywords;
	statement3 %= typeQualifier >> keywords >> '=' >> value;

	Action  action1("statement1", state_);
	Action  action2("statement2", state_);
	Action3 action3("statement3", state_);
	statement %= statement3[action3] | statement2[action2] | statement1[action1];

	IteratorType first = str.begin();
	IteratorType last  = str.end();
	qi::phrase_parse(first, last, statement % ";", AINUR_COMMENTS);

	++first;
	if (first != last && !allEmpty(first, last)) {

		qi::parse(
		    first, last, +('#' >> *(qi::char_ - qi::eol) >> qi::eol | qi::eol | qi::space));

		if (first + 1 != last && !allEmpty(first, last)) {
			IteratorType e = (first + 20 < last) ? first + 20 : last;
			err(AinurState::errLabel(AinurState::ERR_PARSE_FAILED, String(first, e)));
		}
	}

	// Deal with macros
	initMacros();
}

String AinurState::ZERO_CHAR_STRING_(1, ' ');

} // namespace PsimagLite
