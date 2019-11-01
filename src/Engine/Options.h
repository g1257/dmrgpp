#ifndef OPTIONS_H
#define OPTIONS_H
#include "PsimagLite.h"
#include <cctype>

namespace Dmrg {

template<typename InputValidatorType>
class Options {

public:

	typedef typename PsimagLite::String::value_type CharType;
	typedef typename PsimagLite::String::const_iterator StringConstIterator;

	Options(PsimagLite::String label, InputValidatorType& io)
	{
		io.readline(data_, label);
	}

	void operator+=(PsimagLite::String moreData)
	{
		data_ += moreData;
	}

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		ioSerializer.write(label, data_);
	}

	bool isSet(PsimagLite::String what) const
	{
		return isSubstrCaseInsensitive(data_, what);
	}

private:

	template<typename charT>
	struct CaseInsensitiveEqual {

		bool operator()(charT ch1, charT ch2) const
		{
			return std::toupper(ch1) == std::toupper(ch2);
		}
	};

	bool isSubstrCaseInsensitive(const PsimagLite::String& str1,
	                             const PsimagLite::String& str2) const
	{
		StringConstIterator it = std::search(str1.begin(),
		                                     str1.end(),
		                                     str2.begin(),
		                                     str2.end(),
		                                     caseInsensitiveEqual_);
		return (it != str1.end());
	}

	CaseInsensitiveEqual<CharType> caseInsensitiveEqual_;
	PsimagLite::String data_;
};
}
#endif // OPTIONS_H
