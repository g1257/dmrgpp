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
template <typename A, typename ContextType>
typename EnableIf<!TypesEqual<std::vector<std::vector<T> >, A>::True,void>::Type
AinurState::ActionMatrix<T>::operator()(A& attr,
                                        ContextType&,
                                        bool&) const
{
	err("Ainur: Not sure how to handle this\n");
}

template <typename T>
template <typename A, typename ContextType>
typename EnableIf<TypesEqual<std::vector<std::vector<T> >, A>::True,void>::Type
AinurState::ActionMatrix<T>::operator()(A& attr,
                                        ContextType&,
                                        bool&) const
{
	SizeType rows = attr.size();
	if (rows == 0) return;
	SizeType cols = attr[0].size();
	t_.resize(rows, cols);
	for (SizeType i = 0; i < rows; ++i) {
		if (attr[i].size() != cols)
			err("Ainur: Problem reading matrix\n");
		for (SizeType j = 0; j < cols; ++j)
			t_(i, j) = attr[i][j];
	}
}

template <typename T>
template <typename A, typename ContextType>
typename EnableIf<TypesEqual<A,T>::True,void>::Type
AinurState::Action<T>::operator()(A& attr,
                                  ContextType&,
                                  bool&) const
{
	if (name_ == "elipsis") {
		SizeType total = t_.size();
		if (total == 0)
			err("Elipsis cannot be used for this vector, because size is unknown\n");

		for (SizeType i = 0; i < total; ++i)
			t_[i] = attr;
		return;
	}

	t_.push_back(attr);
}

template <typename T>
template <typename A, typename ContextType>
typename EnableIf<!TypesEqual<A,T>::True, void>::Type
AinurState::Action<T>::operator()(A& attr,
                                  ContextType&,
                                  bool&) const
{
	t_ = attr;
}

template <typename A, typename ContextType>
void AinurState::ActionCmplx::operator()(A& attr,
                                         ContextType&,
                                         bool&) const
{
	typedef boost::fusion::vector2<double, double> PairOfDoublesType;
	if (attr.which() == 0) {
		PairOfDoublesType v = boost::get<PairOfDoublesType>(attr);
		t_.push_back(std::complex<double>(
		                 boost::fusion::at_c<0>(v),
		                 boost::fusion::at_c<1>(v)));
	} else {
		t_.push_back(boost::get<double>(attr));
	}
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
	if (!r)
		err("matrix parsing failed\n");

	if (it != value.end())
		std::cerr << "matrix parsing: unmatched part exists\n";
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

	Action<T> actionRows("rows", t);
	Action<T> actionElipsis("elipsis", t);

	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          ruRows [actionRows] | ruElipsis[actionElipsis],
	                          qi::space);

	//check if we have a match
	if (!r)
		err("vector parsing failed\n");

	if (it != value.end())
		std::cerr << "vector parsing: unmatched part exists\n";
}

template<typename T>
void AinurState::convertInternal(std::vector<std::complex<T> >& t,
                                 String value,
                                 typename EnableIf<Loki::TypeTraits<T>::isArith,
                                 int>::Type) const
{ 
	namespace qi = boost::spirit::qi;
#define MY_COMPLEX (qi::double_ >> qi::double_ >>  "i") | qi::double_
	typedef BOOST_TYPEOF(MY_COMPLEX) OneType;
	typedef std::string::iterator IteratorType;

	ActionCmplx actionRows("rows", t);
	IteratorType it = value.begin();


	OneType cmplxN =  MY_COMPLEX;
#define MY_V_OF_COMPLEX "[" >> (cmplxN [actionRows] % ",") >> "]"
	typedef BOOST_TYPEOF(MY_V_OF_COMPLEX) ManyTypes;
	ManyTypes ruRows = MY_V_OF_COMPLEX;
	//	qi::rule<IteratorType, std::complex<T>(), qi::space_type> ruElipsis =
	//	        ruleElipsis<std::complex<T> >();


	//Action<std::complex<T> > actionElipsis("elipsis", t);

	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          ruRows,// | ruElipsis[actionElipsis],
	                          qi::space);

	//check if we have a match
	if (!r)
		err("vector parsing failed\n");

	if (it != value.end())
		std::cerr << "vector parsing: unmatched part exists\n";
}

template void AinurState::convertInternal(Matrix<double>&,String, int) const;

template void AinurState::convertInternal(std::vector<double>&, String, int) const;
template void AinurState::convertInternal(std::vector<SizeType>&, String, int) const;

template void AinurState::convertInternal(std::vector<std::complex<double> >&,
String, int) const;
} // namespace PsimagLite
