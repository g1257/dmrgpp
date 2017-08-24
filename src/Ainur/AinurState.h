#ifndef AINURSTATE_H
#define AINURSTATE_H
#include "../Vector.h"

namespace PsimagLite {

class AinurState {

	typedef Vector<String>::Type VectorStringType;

	enum ErrorEnum
	{
		ERR_PARSE_UNDECLARED,
		ERR_PARSE_DECLARED,
		ERR_READ_UNDECLARED,
		ERR_READ_NO_VALUE
	};

public:

	AinurState()
	{
		assert(ZERO_CHAR_STRING_.length() == 1);
		if (ZERO_CHAR_STRING_[0] != ' ')
			err("Ainur::AinurState should be a singleton\n");

		ZERO_CHAR_STRING_[0] = 0;
	}

	void assign(String k, String v)
	{
		int x = storageIndexByName(k);
		if (x < 0)
			err(errLabel(ERR_PARSE_UNDECLARED, k));

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
			err(errLabel(ERR_READ_UNDECLARED, label));
		assert(static_cast<SizeType>(x) < values_.size());
		String val = values_[x];
		if (isEmptyValue(val))
			err(errLabel(ERR_READ_NO_VALUE, label));

		convertInternal(t, val);
	}

	static bool verbose() { return false; }

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

	void convertInternal(SizeType& t, String label) const
	{
		t = atoi(label.c_str());
	}

	static bool isEmptyValue(String s)
	{
		return (s.length() == 0 || s == ZERO_CHAR_STRING_);
	}

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
		case ERR_READ_UNDECLARED:
			return "FATAL read error: No such label " + key + "\n" +
			        "The label " + key + " must appear in the input file\n";
		case ERR_READ_NO_VALUE:
			return "No value provided for label " + key + "\n" +
			        "The label " + key + " must appear in the input file with " +
			        "a non-empty value\n";
		default:
			return "FATAL Ainur error: Unknown error\n";
		}
	}

	static String ZERO_CHAR_STRING_;
	VectorStringType typesDotified_;
	VectorStringType keys_;
	VectorStringType values_;
};

}
#endif // AINURSTATE_H
