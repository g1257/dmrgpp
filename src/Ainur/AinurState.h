#ifndef AINURSTATE_H
#define AINURSTATE_H
#include "../Vector.h"
#include <cassert>
#include "../PsimagLite.h"
#include <numeric>
#include "../Matrix.h"

namespace PsimagLite {

class AinurState {

	typedef Vector<SizeType>::Type VectorSizeType;
	typedef Vector<String>::Type VectorStringType;
	typedef std::complex<double> ComplexType;

	struct myprint
	{
		template <typename T>
		void operator()(const T &t) const
		{
			std::cout << " --------> " << t << '\n';
		}
	};

	struct ActionCmplx {

		ActionCmplx(String name, std::vector<ComplexType>& t)
		    : name_(name), t_(t)
		{}

		template <typename A, typename ContextType>
		typename EnableIf<TypesEqual<A, double>::True, void>::Type
		operator()(A& attr,
		           ContextType&,
		           bool&) const;

		template <typename A, typename ContextType>
		typename EnableIf<!TypesEqual<A, double>::True, void>::Type
		operator()(A& attr,
		           ContextType&,
		           bool&) const;


		String name_;
		std::vector<ComplexType>& t_;
	};

	template<typename T>
	struct Action {

		Action(String name, std::vector<T>& t)
		    : name_(name), t_(t)
		{}

		template <typename A, typename ContextType>
		typename EnableIf<TypesEqual<A,T>::True, void>::Type
		operator()(A& attr,
		           ContextType&,
		           bool&) const;

		template <typename A, typename ContextType>
		typename EnableIf<!TypesEqual<A,T>::True, void>::Type
		operator()(A& attr,
		           ContextType& context,
		           bool& hit) const;

	private:

		String name_;
		std::vector<T>& t_;
	}; // struct Action

	template<typename T>
	struct ActionMatrix {

		ActionMatrix(String name, Matrix<T>& t)
		    : name_(name), t_(t)
		{}

		template <typename A, typename ContextType>
		typename EnableIf<!TypesEqual<std::vector<std::vector<T> >, A>::True,void>::Type
		operator()(A& attr,
		           ContextType&,
		           bool&) const;

		template <typename A, typename ContextType>
		typename EnableIf<TypesEqual<std::vector<std::vector<T> >, A>::True,void>::Type
		operator()(A& attr,
		           ContextType&,
		           bool&) const;

	private:

		String name_;
		Matrix<T>& t_;
	}; // struct ActionMatrix

public:

	enum ErrorEnum
	{
		ERR_PARSE_UNDECLARED,
		ERR_PARSE_DECLARED,
		ERR_PARSE_FAILED,
		ERR_READ_UNDECLARED,
		ERR_READ_NO_VALUE
	};

	AinurState()
	{
		assert(ZERO_CHAR_STRING_.length() == 1);
		if (ZERO_CHAR_STRING_[0] != ' ')
			err("Ainur::AinurState should be a singleton\n");

		ZERO_CHAR_STRING_[0] = 0;
	}

	void assign(String k, String v);

	void declare(String d, String k, String v)
	{
		assignStorageByName(k);
		SizeType u = 1;
		SizeType last = d.length();
		assert(last > 0);
		if (last > 1 && d[last - 1] == '!') {
			u = 0;
			d = d.substr(0,last - 1);
		}

		typesSpec_.push_back(d);
		values_.push_back(v);
		used_.push_back(u);
	}

	void declare(String d, String k)
	{
		declare(d, k, ZERO_CHAR_STRING_);
	}

	void printUnused(std::ostream& os) const
	{
		SizeType n = used_.size();
		bool flag = false;
		for (SizeType i = 0; i < n; ++i) {
			if (used_[i] > 0) continue;
			flag = true;
			break;
		}

		if (!flag) return;

		os<<"Unused keys:\n";

		if (n != keys_.size())
			err("printUnused: internal error\n");

		for (SizeType i = 0; i < n; ++i) {
			if (used_[i] > 0) continue;
			os<<keys_[i]<<"\n";
		}
	}

	void printAll(std::ostream& os) const
	{
		SizeType n = keys_.size();
		assert(n == values_.size());
		for (SizeType i = 0; i < n; ++i)
			os<<keys_[i]<<" "<<values_[i]<<"\n";
	}

	template<typename SomeType>
	void readValue(SomeType& t, String label) const
	{
		int x = storageIndexByName(label);
		if (x < 0)
			err(errLabel(ERR_READ_UNDECLARED, label));
		assert(static_cast<SizeType>(x) < values_.size());
		String val = values_[x];
		if (isEmptyValue(val))
			err(errLabel(ERR_READ_NO_VALUE, label));

		assert(static_cast<SizeType>(x) < typesSpec_.size());
		convertInternal(t, val);
		used_[x]++;
	}

	static bool verbose() { return false; }

	static String errLabel(ErrorEnum e, String key)
	{
		switch (e) {
		case ERR_PARSE_UNDECLARED:
			return "FATAL parse error: Undeclared " + key + "\n" +
			        "You provided a label in the " +
			        "input file that was not recognized.\n" +
			        "Please check the spelling. If you intended " +
			        "to introduce a temporary label,\nyou must declare " +
			        "it first; please see Ainur input format documentation.\n";
		case ERR_PARSE_DECLARED:
			return "FATAL parse error: Already declared " + key + "\n" +
			        "You tried to re-declare a variable that was already declared.\n" +
			        "If you intended to just provide a value for " + key +
			        " then please remove the declaration word.\n";
		case ERR_PARSE_FAILED:
			return "FATAL parse error: Parsing failed near " + key + "\n" +
			        "This is probably a syntax error.\n";
		case ERR_READ_UNDECLARED:
			return "FATAL read error: No such label " + key + "\n" +
			        "The label " + key + " must appear in the input file\n";
		case ERR_READ_NO_VALUE:
			return "FATAL read error: No value provided for label " + key + "\n" +
			        "The label " + key + " must appear in the input file with " +
			        "a non-empty value\n";
		default:
			return "FATAL Ainur error: Unknown error\n";
		}
	}

private:

	int assignStorageByName(String key)
	{
		int x = storageIndexByName(key);
		if (x >= 0)
			err(errLabel(ERR_PARSE_DECLARED, key));
		keys_.push_back(key);
		return keys_.size() - 1;
	}

	int storageIndexByName(String key) const
	{
		VectorStringType::const_iterator it = std::find(keys_.begin(),
		                                                keys_.end(),
		                                                key);
		if (it == keys_.end())
			return -1;
		return it - keys_.begin();
	}


	template<typename T>
	void convertInternal(T& t,
	                     String value,
	                     typename EnableIf<!Loki::TypeTraits<T>::isIntegral
	                     && !Loki::TypeTraits<T>::isFloat,
	                     int>::Type = 0) const
	{
		err("convertInternal generic type value = " + value + "\n");
	}

	template<typename T>
	void convertInternal(std::vector<T>& t,
	                     String value,
	                     typename EnableIf<Loki::TypeTraits<T>::isArith,
	                     int>::Type = 0) const;

	template<typename T>
	void convertInternal(std::vector<std::complex<T> >& t,
	                     String value,
	                     typename EnableIf<Loki::TypeTraits<T>::isArith,
	                     int>::Type = 0) const;

	template<typename T>
	void convertInternal(Matrix<T>& t,
	                     String value,
	                     typename EnableIf<Loki::TypeTraits<T>::isArith,
	                     int>::Type = 0) const;

	template<typename T>
	void convertInternal(T& t,
	                     String label,
	                     typename EnableIf<Loki::TypeTraits<T>::isIntegral,
	                     int>::Type = 0) const
	{
		t = atoi(label.c_str());
	}

	template<typename T>
	void convertInternal(T& t,
	                     String label,
	                     typename EnableIf<Loki::TypeTraits<T>::isFloat,
	                     int>::Type = 0) const
	{
		t = atof(label.c_str());
	}

	void convertInternal(String& t, String label) const
	{
		SizeType l = label.size();
		if (l > 1 && label[0] == '"' && label[l - 1] == '"') {
			t = (l == 2) ? "" : label.substr(1,l - 2);
			return;
		}

		t = label;
	}

	static bool isEmptyValue(String s)
	{
		return (s.length() == 0 || s == ZERO_CHAR_STRING_);
	}

	static String ZERO_CHAR_STRING_;
	VectorStringType typesSpec_;
	VectorStringType keys_;
	VectorStringType values_;
	mutable VectorSizeType used_;
};

}
#endif // AINURSTATE_H
