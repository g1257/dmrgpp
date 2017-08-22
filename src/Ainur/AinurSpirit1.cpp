#include "AinurSpirit.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace PsimagLite {

template <typename A, typename ContextType>
void Ainur::Action::operator()(A& attr,
                               ContextType&,
                               bool&) const
{
	std::cerr <<" ****************** Action name = "<<name_<<"\n";
	std::cerr << "typeid(A).name() = " << typeid(A).name() << "\n";
	std::cerr << "typeid(ContextType).name() = " << typeid(ContextType).name() << "\n";
	//std::cout << "typeid(Locals).name()     = " << typeid(Locals).name() << "\n";

	std::cerr << "attributes: \n";
	boost::fusion::for_each(attr, myprint());

	std::cerr << "attr = "<<attr<<"\n";
	//std::cout << "hit = "<<hit<<"\n";

	if (name_ == "statement2") {
		// String v1 = boost::fusion::at_c<0>(attr);
		//		String v2 = boost::fusion::at_c<2>(attr);

		// keys_.push_back(v1);
		//		values_.push_back(v2);
	}

}

Ainur::Ainur(String str)
    : dummy_("")
{
	std::cerr<<str<<"\n\n";
#define AINUR_COMMENTS ('#' >> *(qi::char_ - qi::eol) >> qi::eol) | qi::eol | qi::space
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	typedef boost::fusion::vector<std::string, std::string> AttribType;
	typedef boost::fusion::vector<std::string, std::string, std::string> Attrib3Type;
	typedef  boost::spirit::qi::locals<State> LocalsType;

	typedef BOOST_TYPEOF(AINUR_COMMENTS) SkipperType;

	qi::rule<IteratorType, std::string(), qi::unused_type> aToZ;
	qi::rule<IteratorType, std::string(), qi::unused_type> zeroToNine;
	qi::rule<IteratorType, std::string(), qi::unused_type> keywords;
	qi::rule<IteratorType, std::string(), qi::unused_type> value;
	qi::rule<IteratorType, std::string(), qi::unused_type> typeQualifier;
	qi::rule<IteratorType, AttribType, LocalsType, SkipperType> statement1;
	qi::rule<IteratorType, AttribType, LocalsType, SkipperType> statement2;
	qi::rule<IteratorType, Attrib3Type, LocalsType, SkipperType> statement3;
	qi::rule<IteratorType, SkipperType> statement;
	value %= +(qi::char_ - (qi::char_(";") | qi::space | qi::eol));
	aToZ = ascii::char_("a","z") | ascii::char_("A", "Z");
	zeroToNine = ascii::char_("0","9");
	typeQualifier %= +(aToZ | ascii::char_("."));
	keywords = +aToZ >> *(aToZ | zeroToNine | ":");
	statement1   %= keywords >> '=' >> value;
	statement2 %= typeQualifier >> keywords;
	statement3 %= typeQualifier >>  keywords >> '=' >> value;
	statement %= statement3 [Action("statement3")] |
	        statement2 [Action("statement2")] |
	        statement1 [Action("statement1")];

	IteratorType first = str.begin();
	IteratorType last = str.end();
	bool r = qi::phrase_parse(first,
	                          last,
	                          statement % ";",
	                          AINUR_COMMENTS);

	bool finished = (first != last);// fail if we did not get a full match
	std::cout<<"finished="<<finished<<" r= "<<r<<"\n";
}

} // namespace PsimagLite

