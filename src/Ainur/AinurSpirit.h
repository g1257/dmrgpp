#ifndef _AINUR_SPIRIT_H_
#define _AINUR_SPIRIT_H_
#include "../Vector.h"
#include "../TypeToString.h"
#include "../PsimagLite.h"
#include <iostream>
#include <string>
#include <fstream>
#include "AinurState.h"

namespace PsimagLite {

class Ainur {

	struct Action {

		Action(String name, AinurState& state)
		    : name_(name), state_(state)
		{}

		template <typename A, typename ContextType>
		void operator()(A& attr,
		                ContextType& context,
		                bool& hit) const;

	private:

		String name_;
		AinurState& state_;
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
	typedef Vector<char>::Type VectorCharType;

	Ainur(String str);

	String& prefix() { return dummy_; }

	const String& prefix() const { return dummy_; }

	void printUnused(std::ostream& os) const
	{
		os<<"PRINT UNUSED\n";
	}

	void printAll(std::ostream& os) const
	{
		state_.printAll(os);
	}

	template<typename SomeType>
	void readValue(SomeType& t, String label) const
	{
		state_.readValue(t, dummy_ + label);
	}

private:

	String dummy_;
	AinurState state_;
}; // class AinurSpirit

}
#endif // _AINUR_SPIRIT_H_
