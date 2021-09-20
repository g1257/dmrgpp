#include "AinurState.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include "AinurDoubleOrFloat.h"

namespace PsimagLite {

struct MyProxyFor {

	static void convert(long unsigned int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(unsigned int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(long int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(int& t, std::string str)
	{
		t = PsimagLite::atoi(str);
	}

	static void convert(double& t, std::string str)
	{
		t = PsimagLite::atof(str);
	}

	static void convert(float& t, std::string str)
	{
		t = PsimagLite::atof(str);
	}

	template<typename T>
	static void convert(std::complex<T>& t, std::string str)
	{
		t = toComplex<T>(str);
	}

	static void convert(String& t, std::string str)
	{
		t = str;
	}

	template<typename T>
	static void convert(T& t, std::string str)
	{
		throw RuntimeError("complex not implemented yet\n");
	}

private:

	template<typename RealType>
	static std::complex<RealType> toComplex(std::string str)
	{
		typedef std::complex<RealType> ComplexType;
		String buffer;
		bool flag = false;
		const SizeType n = str.length();
		RealType real1 = 0;
		for (SizeType i = 0; i < n; ++i) {
			bool isSqrtMinus1 = (str[i] == 'i');
			if (isSqrtMinus1 && flag)
				throw RuntimeError("Error parsing number " + str + "\n");

			if (isSqrtMinus1) {
				flag = true;
				real1 = atof(buffer.c_str());
				buffer = "";
				continue;
			}

			buffer += str[i];
		}

		return (flag) ? ComplexType(real1, atof(buffer.c_str())) :
		                ComplexType(atof(buffer.c_str()), 0);
	}

};

//---------

boost::spirit::qi::rule<std::string::iterator,
std::vector<std::string>(),
boost::spirit::qi::space_type>
ruleRows()
{
	boost::spirit::qi::rule<std::string::iterator,
	        std::vector<std::string>(),
	        boost::spirit::qi::space_type> myrule =  "[" >> (+~boost::spirit::qi::char_(",[]"))
	                                                        % ',' >> "]";
	return myrule;
}

//---------

void AinurState::assign(String k, String v)
{
	int x = storageIndexByName(k);
	if (x < 0)
		err(errLabel(ERR_PARSE_UNDECLARED, k));

	assert(static_cast<SizeType>(x) < values_.size());

	//if (values_[x] != "")
	//	std::cerr<<"Overwriting label "<<k<<" with "<<v<<"\n";

	values_[x] = v;
}

//---------

template <typename T>
template <typename A, typename ContextType>
void AinurState::ActionMatrix<T>::operator()(A& attr,
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
			MyProxyFor::convert(t_(i, j), attr[i][j]);
	}
}

//---------

template <typename T>
template <typename A, typename ContextType>
void AinurState::Action<T>::operator()(A& attr,
                                       ContextType&,
                                       bool&) const
{
	const SizeType n = attr.size();
	if (n == 2 && attr[1] == "...") {
		const SizeType m = t_.size();
		if (m == 0)
			err("Cannot use ellipsis for vector of unknown size\n");
		MyProxyFor::convert(t_[0], attr[0]);
		for (SizeType i = 1; i < m; ++i)
			t_[i] = t_[0];
		return;
	}

	if (n == 2 && attr[1].length() > 4 && attr[1].substr(0, 4) == "...x") {
		const SizeType m = t_.size();
		const SizeType l = attr[1].length();
		const SizeType mm = PsimagLite::atoi(attr[1].substr(4, l - 4));
		if (m != 0)
			std::cout<<"Resizing vector to "<<mm<<"\n";
		t_.resize(mm);
		MyProxyFor::convert(t_[0], attr[0]);
		for (SizeType i = 1; i < mm; ++i)
			t_[i] = t_[0];
		return;
	}

	t_.resize(n);
	for (SizeType i = 0; i < n; ++i)
		MyProxyFor::convert(t_[i], attr[i]);
}

//---------

template<typename T>
void AinurState::convertInternal(Matrix<T>& t,
                                 String value) const
{
	namespace qi = boost::spirit::qi;
	typedef std::string::iterator IteratorType;
	typedef std::vector<std::string> VectorStringType;
	typedef std::vector<VectorStringType> VectorVectorVectorType;

	IteratorType it = value.begin();
	qi::rule<IteratorType, VectorStringType(), qi::space_type> ruRows = ruleRows();

	qi::rule<IteratorType, VectorVectorVectorType(), qi::space_type> full =
	        "[" >> -(ruRows  % ",") >> "]";

	ActionMatrix<T> actionMatrix("matrix", t);
	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          full [actionMatrix],
	                          qi::space);

	//check if we have a match
	if (!r) {
		err("matrix parsing failed near " +
		    stringContext(it, value.begin(), value.end())
		    + "\n");
	}

	if (it != value.end())
		std::cerr << "matrix parsing: unmatched part exists\n";
}

template<typename T>
void AinurState::convertInternal(std::vector<T>& t,
                                 String value,
                                 typename EnableIf<Loki::TypeTraits<T>::isArith ||
                                 IsComplexNumber<T>::True ||
                                 TypesEqual<T, String>::True,
                                 int>::Type) const
{
	namespace qi = boost::spirit::qi;
	typedef std::string::iterator IteratorType;
	typedef std::vector<std::string> VectorStringType;

	IteratorType it = value.begin();
	qi::rule<IteratorType, VectorStringType(), qi::space_type> ruRows = ruleRows();

	Action<T> actionRows("rows", t);

	bool r = qi::phrase_parse(it,
	                          value.end(),
	                          ruRows [actionRows],
	                          qi::space);

	//check if we have a match
	if (!r)
		err("vector parsing failed near " + stringContext(it, value.begin(), value.end()) + "\n");

	if (it != value.end()) {
		std::cerr << "vector parsing: unmatched part exists near ";
		std::cerr << stringContext(it, value.begin(), value.end())<<"\n";
	}
}

//---------

template void AinurState::convertInternal(Matrix<DoubleOrFloatType>&,String) const;

template void AinurState::convertInternal(Matrix<std::complex<DoubleOrFloatType> >&,
String) const;

template void AinurState::convertInternal(Matrix<String>&, String) const;

template void AinurState::convertInternal(std::vector<DoubleOrFloatType>&, String, int) const;

template void AinurState::convertInternal(std::vector<std::complex<DoubleOrFloatType> >&,
String,
int) const;

template void AinurState::convertInternal(std::vector<SizeType>&, String, int) const;

template void AinurState::convertInternal(std::vector<int>&, String, int) const;

template void AinurState::convertInternal(std::vector<String>&, String, int) const;
} // namespace PsimagLite
