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
	//boost::fusion::for_each(context.attributes, myprint());

	std::cout << "locals: \n";
	//boost::fusion::for_each(context.locals, myprint());

	std::cout << "attr = "<<attr<<"\n";
	std::cout << "hit = "<<hit<<"\n";
}

Ainur::Ainur(String str)
    : dummy_("")
{
	namespace qi = boost::spirit::qi;
	namespace ascii = boost::spirit::ascii;
	typedef boost::fusion::vector<std::string, std::string> AttribType;

	qi::rule<IteratorType, std::string(), qi::space_type> keywords_;
	qi::rule<IteratorType, std::string(), qi::space_type> quoted2_;
	qi::rule<IteratorType, std::string(), qi::space_type> quotedString_;
	qi::rule<IteratorType, AttribType, qi::space_type> statement1_;

	quotedString_ %= qi::lexeme['"' >> +(qi::char_ - '"')  >> '"'];
	quoted2_ %= quotedString_ [Action("lexeme2")];
	keywords_ %= (+(ascii::char_("a","z") || ascii::char_("A", "Z") || ascii::char_('_')));
	statement1_   %= keywords_  [ Action("keywords") ] >> '='
	                                                >> quoted2_ [ Action("quotedString") ];
	Iterator first = str.begin();
	Iterator last = str.end();
	bool r = qi::phrase_parse(first,
	                          last,
	                          statement1_ [Action("statement1")] % ";",
	        qi::space);

	bool finished = (first != last);// fail if we did not get a full match
	std::cout<<"finished="<<finished<<" r= "<<r<<"\n";
}

} // namespace PsimagLite

