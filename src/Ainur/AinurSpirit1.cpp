#include "AinurSpirit.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace PsimagLite {

template <typename A, typename ContextType>
void Ainur::Action::operator()(A& attr, const ContextType&, bool&) const
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
#define AINUR_COMMENTS ('#' >> *(qi::char_ - qi::eol) >> qi::eol) | qi::eol
	//| qi::space
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	typedef boost::fusion::vector<std::string, std::string> AttribType;
	typedef boost::fusion::vector<std::string, std::string, std::string> Attrib2Type;
	typedef BOOST_TYPEOF(AINUR_COMMENTS) SkipperType;

	qi::rule<IteratorType, std::string(), qi::unused_type> keywords;
	qi::rule<IteratorType, std::string(), qi::unused_type> value;
	//qi::rule<IteratorType, std::string()> quotedString;
	qi::rule<IteratorType, AttribType, SkipperType> statement1;
	qi::rule<IteratorType, Attrib2Type, SkipperType> statement2;
	qi::rule<IteratorType, SkipperType> statement;
	//quotedString %= qi::lexeme['"' >> +(qi::char_ - '"')  >> '"'];
	value %= +(qi::char_ - (qi::char_(";") | qi::space | qi::eol));
	keywords %= (+(ascii::char_("a","z") |
	               ascii::char_("A", "Z") |
	               ascii::char_(".") |
	               ascii::char_(":")));
	statement1   %= keywords >> '=' > value;
	statement2 %= keywords >> +(qi::space) >> value;
	statement %= (statement1 [Action("statement1")] | statement2 [Action("statement2")]);

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

