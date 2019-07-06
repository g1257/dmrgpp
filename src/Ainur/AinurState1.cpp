#include "AinurState.h"
#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include "AinurDoubleOrFloat.h"

namespace PsimagLite {

template<typename T>
boost::spirit::qi::rule<std::string::iterator,
std::vector<T>(),
boost::spirit::qi::space_type>
ruleRows();

template<>
boost::spirit::qi::rule<std::string::iterator,
std::vector<DoubleOrFloatType>(),
boost::spirit::qi::space_type>
ruleRows<DoubleOrFloatType>()
{
	return "[" >> -(boost::spirit::DoubleOrFloatUnderscore % ",") >> "]";
}

template<typename SomeRealType>
class MyClass {

public:

	typedef std::complex<SomeRealType> value_type;

	MyClass() : t1_(0), t2_(0), counter_(0) {}

	explicit MyClass(const SomeRealType& r) : t1_(r), t2_(0), counter_(0) {}

	value_type toComplex() const
	{
		return (counter_ < 2) ? t1_ : value_type(t1_, t2_);
	}

	int end() const { return 0; }

	void insert(int, const SomeRealType& val)
	{
		if (counter_ == 0)
			t1_ = val;
		else if (counter_ == 1)
			t2_ = val;
		if (counter_ > 1)
			err("Wrong complex\n");
		++counter_;
	}

	void insert(int, const value_type& val)
	{
		err("Wrong complex ...\n");
	}

private:

	SomeRealType t1_;
	SomeRealType t2_;
	SizeType counter_;
};

template<typename ComplexOrRealType, bool>
struct MyProxyFor {
	typedef ComplexOrRealType Type;

	static void copy(ComplexOrRealType& dest, const ComplexOrRealType& src)
	{
		dest = src;
	}

	static void copy(std::vector<ComplexOrRealType>& dest,
	                 const std::vector<ComplexOrRealType>& src)
	{
		dest = src;
	}
};

template<typename ComplexOrRealType>
struct MyProxyFor<ComplexOrRealType, true> {

	typedef MyClass<typename Real<ComplexOrRealType>::Type> Type;

	static void copy(std::vector<ComplexOrRealType>& dest,
	                 const std::vector<Type>& src)
	{
		dest.clear();
		const SizeType n = src.size();
		if (n == 0) return;
		dest.resize(n);
		for (SizeType i = 0; i < n; ++i)
			dest[i] = src[i].toComplex();
	}

	static void copy(std::vector<ComplexOrRealType>& dest,
	                 const Type& src)
	{
		dest.push_back(src.toComplex());
	}
};

template<>
boost::spirit::qi::rule<std::string::iterator,
std::vector<MyClass<DoubleOrFloatType> >(),
boost::spirit::qi::space_type>
ruleRows<MyClass<DoubleOrFloatType> >()
{
	auto fl = boost::spirit::DoubleOrFloatUnderscore;
	boost::spirit::qi::rule<std::string::iterator,
	std::vector<MyClass<DoubleOrFloatType> >(),
	boost::spirit::qi::space_type> myrule =  "[" >> -( ( fl  | "(" >> fl >> "," >> fl >> ")" )
	            % ",") >> "]";
	return myrule;
}

template<>
boost::spirit::qi::rule<std::string::iterator,
std::vector<SizeType>(),
boost::spirit::qi::space_type>
ruleRows<SizeType>()
{
	return "[" >> -(boost::spirit::int_ % ",") >> "]";
}

template<>
boost::spirit::qi::rule<std::string::iterator,
std::vector<int>(),
boost::spirit::qi::space_type>
ruleRows<int>()
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
DoubleOrFloatType(),
boost::spirit::qi::space_type>
ruleElipsis<DoubleOrFloatType>()
{
	return  "[" >> boost::spirit::DoubleOrFloatUnderscore  >> "," >> "..." >> "]";
}

template<>
boost::spirit::qi::rule<std::string::iterator,
MyClass<DoubleOrFloatType>(),
boost::spirit::qi::space_type>
ruleElipsis<MyClass<DoubleOrFloatType> >()
{
	auto fl = boost::spirit::DoubleOrFloatUnderscore;
	boost::spirit::qi::rule<std::string::iterator,
	MyClass<DoubleOrFloatType>(),
	boost::spirit::qi::space_type> myrule =  "[" >> (fl |
	                (fl >> "i" >> fl) |
	                ("i" >> fl))
	            >> "," >> "..." >> "]";
	return myrule;
}

template<>
boost::spirit::qi::rule<std::string::iterator,
SizeType(),
boost::spirit::qi::space_type>
ruleElipsis<SizeType>()
{
	return "[" >> boost::spirit::int_ >> "," >> "..." >> "]";
}

template<>
boost::spirit::qi::rule<std::string::iterator,
int(),
boost::spirit::qi::space_type>
ruleElipsis<int>()
{
	return "[" >> boost::spirit::int_ >> "," >> "..." >> "]";
}

//---------

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
	MyProxyFor<T, IsComplexNumber<T>::True>::copy(t_, attr);
}

template <typename A, typename ContextType>
typename EnableIf<TypesEqual<A, DoubleOrFloatType>::True, void>::Type
AinurState::ActionCmplx::operator()(A& attr,
                                    ContextType&,
                                    bool&) const
{
	if (name_ == "rows2")
		t_.push_back(std::complex<DoubleOrFloatType>(0.0, attr));
	else if (name_ == "rows3")
		t_.push_back(attr);
	else
		err("Error parsing complex number\n");
}

template <typename A, typename ContextType>
typename EnableIf<!TypesEqual<A, DoubleOrFloatType>::True, void>::Type
AinurState::ActionCmplx::operator()(A& attr,
                                    ContextType&,
                                    bool&) const
{
	if (name_ == "rows1") {
		t_.push_back(std::complex<DoubleOrFloatType>(
		                 boost::fusion::at_c<0>(attr),
		                 boost::fusion::at_c<1>(attr)));
	} else {
		err("Error parsing complex number\n");
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
                                 IsComplexNumber<T>::True,
                                 int>::Type) const
{
	namespace qi = boost::spirit::qi;
	typedef std::string::iterator IteratorType;
	typedef typename MyProxyFor<T, IsComplexNumber<T>::True>::Type MyProxyForType;
	typedef std::vector<MyProxyForType> LocalVectorType;

	IteratorType it = value.begin();
	qi::rule<IteratorType, LocalVectorType(), qi::space_type> ruRows = ruleRows<MyProxyForType>();
	qi::rule<IteratorType, MyProxyForType(), qi::space_type> ruElipsis =
	        ruleElipsis<MyProxyForType>();

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

template void AinurState::convertInternal(Matrix<DoubleOrFloatType>&,String, int) const;

template void AinurState::convertInternal(std::vector<DoubleOrFloatType>&, String, int) const;

template void AinurState::convertInternal(std::vector<SizeType>&, String, int) const;

template void AinurState::convertInternal(std::vector<int>&, String, int) const;

template void AinurState::convertInternal(std::vector<std::complex<DoubleOrFloatType> >&,
String, int) const;
} // namespace PsimagLite
