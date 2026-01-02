#ifndef CMPLX_OR_REAL_H
#define CMPLX_OR_REAL_H
#include "Vector.h"

namespace PsimagLite {

template <typename RealType, int> class CpmlxOrReal { };

template <typename RealType> class CpmlxOrReal<RealType, 1> {

public:

	typedef std::complex<RealType> ComplexType;

	CpmlxOrReal(String t)
	    : isImag_(false)
	    , value_(0)
	{

		const SizeType len = t.length();
		String         buffer;
		for (SizeType i = 0; i < len; ++i) {
			if (t[i] == 'i') {
				if (isImag_)
					throw RuntimeError("More than one i found\n");
				isImag_ = true;
				continue;
			}

			buffer += t[i];
		}

		if (isImag_ && buffer == "")
			value_ = 1;
		else
			value_ = PsimagLite::atof(buffer);
	}

	ComplexType value() const { return (isImag_) ? ComplexType(0, value_) : value_; }

private:

	bool     isImag_;
	RealType value_;
};

template <typename RealType> class CpmlxOrReal<RealType, 0> {

public:

	CpmlxOrReal(String t)
	    : value_(0)
	{
		bool isComplex = (t.find("i") != String::npos);
		if (isComplex)
			throw RuntimeError("i \\equiv sqrt(-1) found but code path is real\n");
		value_ = PsimagLite::atof(t);
	}

	RealType value() const { return value_; }

private:

	RealType value_;
};

} // namespace PsimagLite
#endif // CMPLX_OR_REAL_H
