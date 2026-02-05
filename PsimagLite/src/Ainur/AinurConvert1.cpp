#include "AinurComplex.hh"
#include "AinurConvert.hh"
#include "AinurDoubleOrFloat.h"
#include <boost/config/warning_disable.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <boost/spirit/include/qi.hpp>

namespace PsimagLite {

boost::spirit::qi::
    rule<std::string::iterator, std::vector<std::string>(), boost::spirit::qi::space_type>
    ruleRows()
{
	boost::spirit::qi::
	    rule<std::string::iterator, std::vector<std::string>(), boost::spirit::qi::space_type>
	        myrule = "[" >> (+~boost::spirit::qi::char_(",[]")) % ',' >> "]";
	return myrule;
}

template <typename T>
template <typename A, typename ContextType>
void AinurConvert::ActionMatrix<T>::operator()(A& attr, ContextType&, bool&) const
{
	SizeType rows = attr.size();
	if (rows == 0)
		return;
	SizeType cols = attr[0].size();
	t_.resize(rows, cols);
	for (SizeType i = 0; i < rows; ++i) {
		if (attr[i].size() != cols)
			err("Ainur: Problem reading matrix\n");
		for (SizeType j = 0; j < cols; ++j) {
			String attrAfter = ainurMacros_.valueFromFunction(attr[i][j]);
			AinurComplex::convert(t_(i, j), attrAfter);
		}
	}
}

//---------

template <typename T>
template <typename A, typename ContextType>
void AinurConvert::Action<T>::operator()(A& attr, ContextType&, bool&) const
{
	const SizeType n = attr.size();
	if (n == 2 && attr[1] == "...") {
		const SizeType m = t_.size();
		if (m == 0)
			err("Cannot use ellipsis for vector of unknown size\n");
		AinurComplex::convert(t_[0], attr[0]);
		for (SizeType i = 1; i < m; ++i)
			t_[i] = t_[0];
		return;
	}

	if (n == 2 && attr[1].length() > 4 && attr[1].substr(0, 4) == "...x") {
		const SizeType m  = t_.size();
		const SizeType l  = attr[1].length();
		const SizeType mm = PsimagLite::atoi(attr[1].substr(4, l - 4));
		if (m != 0)
			std::cout << "Resizing vector to " << mm << "\n";
		t_.resize(mm);
		AinurComplex::convert(t_[0], attr[0]);
		for (SizeType i = 1; i < mm; ++i)
			t_[i] = t_[0];
		return;
	}

	t_.resize(n);
	for (SizeType i = 0; i < n; ++i) {
		String val = ainurMacros_.valueFromFunction(attr[i]);
		AinurComplex::convert(t_[i], val);
	}
}

//---------

template <typename T> void AinurConvert::convert(Matrix<T>& t, const AinurVariable& ainurVariable)
{
	namespace qi                 = boost::spirit::qi;
	using IteratorType           = std::string::iterator;
	using VectorStringType       = std::vector<std::string>;
	using VectorVectorVectorType = std::vector<VectorStringType>;

	String value = ainurVariable.value;

	IteratorType                                               it     = value.begin();
	qi::rule<IteratorType, VectorStringType(), qi::space_type> ruRows = ruleRows();

	qi::rule<IteratorType, VectorVectorVectorType(), qi::space_type> full
	    = "[" >> -(ruRows % ",") >> "]";

	ActionMatrix<T> actionMatrix("matrix", t, ainurMacros_);
	bool            r = qi::phrase_parse(it, value.end(), full[actionMatrix], qi::space);

	// check if we have a match
	if (!r) {
		err("matrix parsing failed near " + stringContext(it, value.begin(), value.end())
		    + "\n");
	}

	if (it != value.end())
		std::cerr << "matrix parsing: unmatched part exists\n";
}

template <typename T>
void AinurConvert::convert(
    std::vector<T>&      t,
    const AinurVariable& ainurVariable,
    typename EnableIf<Loki::TypeTraits<T>::isArith || IsComplexNumber<T>::True
                          || TypesEqual<T, String>::True,
                      int>::Type)
{
	namespace qi           = boost::spirit::qi;
	using IteratorType     = std::string::iterator;
	using VectorStringType = std::vector<std::string>;

	String                                                     value  = ainurVariable.value;
	IteratorType                                               it     = value.begin();
	qi::rule<IteratorType, VectorStringType(), qi::space_type> ruRows = ruleRows();

	Action<T> actionRows("rows", t, ainurMacros_);

	bool r = qi::phrase_parse(it, value.end(), ruRows[actionRows], qi::space);

	// check if we have a match
	if (!r)
		err("vector parsing failed near " + stringContext(it, value.begin(), value.end())
		    + "\n");

	if (it != value.end()) {
		std::cerr << "vector parsing: unmatched part exists near ";
		std::cerr << stringContext(it, value.begin(), value.end()) << "\n";
	}
}

//---------

template void AinurConvert::convert(Matrix<DoubleOrFloatType>&, const AinurVariable&);

template void AinurConvert::convert(Matrix<std::complex<DoubleOrFloatType>>&, const AinurVariable&);

template void AinurConvert::convert(Matrix<String>&, const AinurVariable&);

template void AinurConvert::convert(std::vector<DoubleOrFloatType>&, const AinurVariable&, int);

template void
AinurConvert::convert(std::vector<std::complex<DoubleOrFloatType>>&, const AinurVariable&, int);

template void AinurConvert::convert(std::vector<SizeType>&, const AinurVariable&, int);

template void AinurConvert::convert(std::vector<int>&, const AinurVariable&, int);

template void AinurConvert::convert(std::vector<String>&, const AinurVariable&, int);
} // namespace PsimagLite
