#ifndef AINURSTATE_H
#define AINURSTATE_H
#include "../Vector.h"
#include <cassert>
#include "../PsimagLite.h"
#include <numeric>
#include "../Matrix.h"
#include "AinurDoubleOrFloat.h"
#include "AinurConvert.hh"

namespace PsimagLite {

class AinurState {

public:

	typedef Vector<SizeType>::Type VectorSizeType;
	typedef Vector<String>::Type VectorStringType;
	typedef std::complex<DoubleOrFloatType> ComplexType;

	struct myprint
	{
		template <typename T>
		void operator()(const T &t) const
		{
			std::cout << " --------> " << t << '\n';
		}
	};


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
		//		if (ZERO_CHAR_STRING_[0] != ' ')
		//			err("Ainur::AinurState should be a singleton\n");

		ZERO_CHAR_STRING_[0] = 0;
	}

	void assign(String k, String v)
	{
		int x = storageIndexByName(k);
		if (x < 0)
			err(errLabel(ERR_PARSE_UNDECLARED, k));

		assert(static_cast<SizeType>(x) < values_.size());

		//if (values_[x] != "")
		//	std::cerr<<"Overwriting label "<<k<<" with "<<v<<"\n";

		values_[x] = v;
	}

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

	void initMacros()
	{
		expandMacrosRecursively();
		installNativeMacros();
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
		AinurConvert::convertInternal(t, val);
		used_[x]++;
	}

	template<typename SomeMapType>
	void setMap(SomeMapType& map) const
	{
		const SizeType n = keys_.size();
		assert(n == values_.size());
		for (SizeType i = 0; i < n; ++i)
			if (used_[i]) map[keys_[i]] = values_[i];
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



	static bool isEmptyValue(String s)
	{
		return (s.length() == 0 || s == ZERO_CHAR_STRING_);
	}

	void installNativeMacros()
	{
		declare("function", "AinurFromFile", "AinurFromFile");
	}

	void expandMacrosRecursively()
	{
		static const SizeType avoidInfMax = 100;
		SizeType avoidInfCounter = 0;
		while (expandMacros()) {
			if (avoidInfCounter++ > avoidInfMax) {
				err("Recursion limit of " + ttos(avoidInfMax) + " exceeded.\n");
			}
		}
	}

	bool expandMacros()
	{
		const SizeType n = keys_.size();
		assert(n == values_.size());
		bool atLeastOneValueHasMacro = false;
		for (SizeType i = 0; i < n; ++i) {
			if (!used_[i]) continue;
			std::pair<bool, PsimagLite::String> macro = expandOneValue(values_[i]);

			if (!macro.first) continue;

			values_[i] = macro.second;
			atLeastOneValueHasMacro = true;
		}

		return atLeastOneValueHasMacro;
	}

	// \[a-zA-Z]+
	std::pair<bool, PsimagLite::String> expandOneValue(PsimagLite::String value) const
	{
		const SizeType n = value.length();
		PsimagLite::String macroName;
		PsimagLite::String retString;
		bool hasAtLeastOneMacro = false;
		SizeType status = 0; // 0 = outsite a macro, 1 = inside a macro
		for (SizeType i = 0; i < n; ++i) {
			const char c = value[i];
			if (c == '\\') {
				status = 1;
				continue;
			}

			if (status == 0) {
				retString += c;
				continue;
			}

			if (isValidCharForMacroName(c)) {
				macroName += c;
			} else {

				int x = storageIndexByName(macroName);
				if (x < 0) err("No macro named " + macroName + "\n");

				assert(static_cast<SizeType>(x) < values_.size());
				retString += unquote(values_[x]) + c;

				macroName = "";
				hasAtLeastOneMacro = true;
				status = 0;
			}
		}

		return std::pair<bool, PsimagLite::String>(hasAtLeastOneMacro, unquote(retString));
	}

	static bool isValidCharForMacroName(char c)
	{
		const bool b1 = (c > 96 && c < 123);
		const bool b2 = (c > 64 && c < 91);
		const bool b3 = (c > 47 && c < 58);

		return (b1 || b2 || b3 || c == '_');
	}

	static PsimagLite::String unquote(PsimagLite::String str)
	{
		if (str.length() == 0) return str;
		if (str[0] == '"') str = str.substr(1, str.length() - 1);
		if (str.length() == 0) return str;
		if (str[str.length() - 1] == '"') str = str.substr(0, str.length() - 1);
		return str;
	}

	static String ZERO_CHAR_STRING_;
	VectorStringType typesSpec_;
	VectorStringType keys_;
	VectorStringType values_;
	mutable VectorSizeType used_;
};

}
#endif // AINURSTATE_H
