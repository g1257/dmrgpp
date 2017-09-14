#include "AinurState.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

namespace PsimagLite {

void AinurState::assign(String k, String v)
{
	int x = storageIndexByName(k);
	if (x < 0)
		err(errLabel(ERR_PARSE_UNDECLARED, k));

	assert(static_cast<SizeType>(x) < values_.size());
	values_[x] = v;
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
	qi::rule<IteratorType, LocalVectorType(), qi::space_type> ruRows =
	        "[" >> -(qi::double_ % ",") >> "]";
	qi::rule<IteratorType, LocalVectorVectorType(), qi::space_type> full =
	        "[" >> -(ruRows  % ",") >> "]";

	LocalVectorVectorType data;
	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          full,
	                          qi::space,
	                          data);

	//check if we have a match
	std::cerr << "matched: " << std::boolalpha << r << '\n';
	if (it != value.end())
		std::cerr << "unmatched part exists\n";

	SizeType rows = data.size();
	if (rows == 0) return;
	SizeType cols = data[0].size();
	t.resize(rows, cols);
	for (SizeType i = 0; i < rows; ++i) {
		if (data[i].size() != cols)
			err("Error in matrix at row " + ttos(i) + "\n");

		for (SizeType j = 0; j < cols; ++j) {
			t(i,j) = data[i][j];
		}
	}
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
	qi::rule<IteratorType, LocalVectorType(), qi::space_type> ruRows =
	        "[" >> -(qi::double_ % ",") >> "]";

	//LocalVectorType data;
	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          ruRows,
	                          qi::space,
	                          t);

	//check if we have a match
	std::cerr << "matched: " << std::boolalpha << r << '\n';
	if (it != value.end())
		std::cerr << "unmatched part exists\n";
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
