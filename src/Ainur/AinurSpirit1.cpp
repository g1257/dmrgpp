#include "AinurSpirit.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace PsimagLite {

template <typename A, typename ContextType, typename PassType>
void Ainur::Action::operator()(A& attr,
                            ContextType& context,
                            PassType hit) const
{
	std::cout <<" ****************** Action name = "<<name_<<"\n";
	std::cout << "typeid(A).name() = " << typeid(A).name() << "\n";
	std::cout << "typeid(ContextType).name() = " << typeid(ContextType).name() << "\n";
	//std::cout << "typeid(Locals).name()     = " << typeid(Locals).name() << "\n";

	std::cout << "attributes: \n";
	// boost::fusion::for_each(context.attributes, myprint());

	std::cout << "locals: \n";
	//boost::fusion::for_each(context.locals, Ainur::myprint());

	std::cout << "attr = "<<attr<<"\n";
	std::cout << "hit = "<<hit<<"\n";
}

Ainur::Ainur(String str)
    : dummy_("")
{
#define AINUR_COMMENTS ('#' >> *(qi::char_ - qi::eol) >> qi::eol) | qi::eol | qi::space
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	typedef boost::fusion::vector<std::string, std::string> AttribType;
	typedef BOOST_TYPEOF(AINUR_COMMENTS) SkipperType;

	qi::rule<IteratorType, std::string()> keywords;
	qi::rule<IteratorType, std::string()> quoted2;
	qi::rule<IteratorType, std::string()> quotedString;
	qi::rule<IteratorType, AttribType, SkipperType> statement1;

	quotedString %= qi::lexeme['"' >> +(qi::char_ - '"')  >> '"'];
	quoted2 %= quotedString [Action("lexeme2")];
	keywords %= (+(ascii::char_("a","z") || ascii::char_("A", "Z")));
	statement1   %= keywords [ Action("keywords") ] >> '='
	                                                >> quoted2 [ Action("quotedString") ];
	IteratorType first = str.begin();
	IteratorType last = str.end();
	bool r = qi::phrase_parse(first,
	                          last,
	                          statement1 [Action("statement1")] % ";",
	        AINUR_COMMENTS);

	bool finished = (first != last);// fail if we did not get a full match
	std::cout<<"finished="<<finished<<" r= "<<r<<"\n";
}

} // namespace PsimagLite

