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

	typedef Vector<String>::Type VectorStringType;

	class State {

	public:

		State()
		{
			assert(ZERO_CHAR_STRING_.length() == 1);
			if (ZERO_CHAR_STRING_[0] != ' ')
				err("Ainur::State should be a singleton\n");

			ZERO_CHAR_STRING_[0] = 0;
		}

		void assign(String k, String v)
		{
			int x = storageIndexByName(k);
			if (x < 0)
				err("Undeclared " + k + " in this scope\n");

			assert(static_cast<SizeType>(x) < values_.size());
			values_[x] = v;
		}

		void declare(String d, String k)
		{
			assignStorageByName(k);
			typesDotified_.push_back(d);
			values_.push_back(ZERO_CHAR_STRING_);
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
				err("No such label " + label + "\n");
			assert(static_cast<SizeType>(x) < values_.size());
			String val = values_[x];
			if (isEmptyValue(val))
				err("No value provided for label " + label + "\n");

			convertInternal(t, val);
		}

		static bool verbose() { return false; }

	private:

		int assignStorageByName(String key)
		{
			int x = storageIndexByName(key);
			if (x >= 0)
				err("Already in scope " + key + "\n");
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

		void convertInternal(SizeType& t, String label) const
		{
			t = atoi(label.c_str());
		}

		static bool isEmptyValue(String s)
		{
			return (s.length() == 0 || s == ZERO_CHAR_STRING_);
		}

		static String ZERO_CHAR_STRING_;
		VectorStringType typesDotified_;
		VectorStringType keys_;
		VectorStringType values_;
	};

	struct Action {

		Action(String name, State& state)
		    : name_(name), state_(state)
		{}

		template <typename A, typename ContextType>
		void operator()(A& attr,
		                ContextType& context,
		                bool& hit) const;

	private:

		String name_;
		State& state_;
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
		state_.readValue(t, label);
	}

private:

	String dummy_;
	State state_;
}; // class AinurSpirit

}
#endif // _AINUR_SPIRIT_H_
