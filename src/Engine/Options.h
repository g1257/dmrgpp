#ifndef OPTIONS_H
#define OPTIONS_H
#include "PsimagLite.h"

namespace Dmrg {

template<typename InputValidatorType>
class Options {

public:

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
		return (data_.find(what) != PsimagLite::String::npos);
	}

private:

	PsimagLite::String data_;
};
}
#endif // OPTIONS_H
