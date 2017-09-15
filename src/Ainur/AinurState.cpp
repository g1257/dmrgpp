#include "AinurState.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace PsimagLite {

template<typename T>
boost::spirit::qi::rule<std::string::iterator,
std::vector<T>(),
boost::spirit::qi::space_type>
ruleRows();

template<>
boost::spirit::qi::rule<std::string::iterator,
std::vector<double>(),
boost::spirit::qi::space_type>
ruleRows<double>()
{
	return "[" >> -(boost::spirit::double_ % ",") >> "]";
}

template<>
boost::spirit::qi::rule<std::string::iterator,
std::vector<SizeType>(),
boost::spirit::qi::space_type>
ruleRows<SizeType>()
{
	return "[" >> -(boost::spirit::int_ % ",") >> "]";
}

//---------
template<typename T>
boost::spirit::qi::rule<std::string::iterator,
T(),
boost::spirit::qi::space_type>
ruleElipsis();

template<>
boost::spirit::qi::rule<std::string::iterator,
double(),
boost::spirit::qi::space_type>
ruleElipsis<double>()
{
	return  "[" >> boost::spirit::double_  >> "," >> "..." >> "]";
}

template<>
boost::spirit::qi::rule<std::string::iterator,
SizeType(),
boost::spirit::qi::space_type>
ruleElipsis<SizeType>()
{
	return "[" >> boost::spirit::int_ >> "," >> "..." >> "]";
}

void AinurState::assign(String k, String v)
{
	int x = storageIndexByName(k);
	if (x < 0)
		err(errLabel(ERR_PARSE_UNDECLARED, k));

	assert(static_cast<SizeType>(x) < values_.size());
	values_[x] = v;
}

template <typename T>
template <typename A, typename ContextType> void
//typename EnableIf<TypesEqual<A,T>::True,void>::Type
AinurState::ActionMatrix<T>::operator()(A& attr,
                                  ContextType&,
                                  bool&) const
{
	err("Action Matrix\n");
}

template <typename T>
template <typename A, typename ContextType>
typename EnableIf<TypesEqual<A,T>::True,void>::Type
AinurState::Action<T>::operator()(A& attr,
                                  ContextType&,
                                  bool&) const
{
	t_.push_back(attr);
}

template <typename T>
template <typename A, typename ContextType>
typename EnableIf<!TypesEqual<A,T>::True,void>::Type
AinurState::Action<T>::operator()(A& attr,
                                  ContextType&,
                                  bool&) const
{
	if (true) {
		std::cerr <<" ****************** Action name = "<<name_<<"\n";
		std::cerr << "typeid(A).name() = " << typeid(A).name() << "\n";
		std::cerr << "typeid(ContextType).name() = " << typeid(ContextType).name() << "\n";
		//std::cout << "typeid(Locals).name()     = " << typeid(Locals).name() << "\n";

		//std::cerr << "attributes: \n";
		//boost::fusion::for_each(attr, myprint());

		// std::cerr << "attr = "<<attr<<"\n";
		//std::cout << "hit = "<<hit<<"\n";
	}

	t_ = attr;

	//	if (name_ == "statement1") {
	//		String v1 = boost::fusion::at_c<0>(attr);
	//		String v2 = boost::fusion::at_c<1>(attr);
	//		state_.assign(v1, v2);
	//	} else if (name_ == "statement2") {
	//		String v1 = boost::fusion::at_c<0>(attr);
	//		String v2 = boost::fusion::at_c<1>(attr);
	//		state_.declare(v1, v2);
	//	}
}

template<typename T>
void AinurState::convertInternal(Matrix<T>& t,
                                 String value,
                                 typename EnableIf<Loki::TypeTraits<T>::isArith,
                                 int>::Type) const
{
	namespace qi = boost::spirit::qi;
	typedef std::string::iterator IteratorType;
	typedef std::vector<T> LocalVectorType;
	typedef std::vector<LocalVectorType> LocalVectorVectorType;

	IteratorType it = value.begin();
	qi::rule<IteratorType, LocalVectorType(), qi::space_type> ruRows = ruleRows<T>();

	qi::rule<IteratorType, LocalVectorVectorType(), qi::space_type> full =
	        "[" >> -(ruRows  % ",") >> "]";

	ActionMatrix<T> actionMatrix("matrix", t);
	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          full [actionMatrix],
	                          qi::space);

	//check if we have a match
	std::cerr << "matrix parsing: matched: " << std::boolalpha << r << '\n';
	if (it != value.end())
		std::cerr << "matrix parsing: unmatched part exists\n";

//	SizeType rows = data.size();
//	if (rows == 0) return;
//	SizeType cols = data[0].size();
//	t.resize(rows, cols);
//	for (SizeType i = 0; i < rows; ++i) {
//		if (data[i].size() != cols)
//			err("Error in matrix at row " + ttos(i) + "\n");

//		for (SizeType j = 0; j < cols; ++j) {
//			t(i,j) = data[i][j];
//		}
//	}
}

template<typename T>
void AinurState::convertInternal(std::vector<T>& t,
                                 String value,
                                 typename EnableIf<Loki::TypeTraits<T>::isArith,
                                 int>::Type) const
{
	namespace qi = boost::spirit::qi;
	typedef std::string::iterator IteratorType;
	typedef std::vector<T> LocalVectorType;

	IteratorType it = value.begin();
	qi::rule<IteratorType, LocalVectorType(), qi::space_type> ruRows = ruleRows<T>();
	qi::rule<IteratorType, T(), qi::space_type> ruElipsis = ruleElipsis<T>();

	//SizeType oldSize = t.size();

	Action<T> actionRows("rows", t);
	Action<T> actionElipsis("elipsis", t);

	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          ruRows [actionRows] | ruElipsis[actionElipsis],
	                          qi::space);

	//check if we have a match
	std::cerr << "vector parsing: matched: " << std::boolalpha << r << '\n';
	if (it != value.end())
		std::cerr << "vector parsing: unmatched part exists\n";
}

template<typename T>
void AinurState::convertInternal(std::vector<std::complex<T> >& t,
                                 String value,
                                 typename EnableIf<Loki::TypeTraits<T>::isArith,
                                 int>::Type) const
{
	err("Ainur: Reading vector of complex is not supported yet (sorry)\n");
}

template void AinurState::convertInternal(Matrix<double>&,String, int) const;

template void AinurState::convertInternal(std::vector<double>&, String, int) const;
template void AinurState::convertInternal(std::vector<SizeType>&, String, int) const;

template void AinurState::convertInternal(std::vector<std::complex<double> >&,
String, int) const;
} // namespace PsimagLite
