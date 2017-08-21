#ifndef _AINUR_SPIRIT_H_
#define _AINUR_SPIRIT_H_
#include "../Vector.h"
#include "../TypeToString.h"
#include "../PsimagLite.h"
#include <iostream>
#include <string>
#include <fstream>

namespace PsimagLite {

class Ainur {

	struct Action {

		Action(String name)
		    : name_(name)
		{}

		template <typename A, typename ContextType> //, typename PassType>
		void operator()(A& attr,
		                const ContextType& context,
		                bool& hit) const;

	private:

		String name_;
	};

	struct myprint
	{
		template <typename T>
		void operator()(const T &t) const
		{
			std::cout << " --------> " << t << '\n';
		}
	};

public:

	typedef std::string::iterator IteratorType;
	typedef Vector<String>::Type VectorStringType;
	typedef Vector<char>::Type VectorCharType;
	//typedef AinurStatements AinurStatementsType;
	//typedef AinurStatementsType::AinurLexicalType AinurLexicalType;

	Ainur(String str);

	String& prefix() { return dummy_; }

	const String& prefix() const { return dummy_; }

	void printUnused(std::ostream& os) const
	{
		os<<"PRINT UNUSED\n";
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
		std::cerr<<"readValue called for label="<<label<<"\n";
		err("Ainur isn't ready, throwing...\n");
	}

private:

	String dummy_;
	VectorStringType keys_;
	VectorStringType values_;
}; // class AinurSpirit

}
#endif // _AINUR_SPIRIT_H_
