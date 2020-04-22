#ifndef PVECTOR_H
#define PVECTOR_H
#include "Vector.h"
#include <cstdlib>

namespace Dmrg {

template<typename VectorWithOffsetType>
class Pvector {

public:

	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	// |P0>=(c?0[0]'*c?0[1]' +  c?1[0]'*c?1[1] - c?0[1]'*c?0[0] - c?1[1]'*c?1[0])|gs>*weight
	// The weight is optional
	Pvector(PsimagLite::String str)
	    : str_(str), weight_(1.0)
	{
		// find the weight first
		SizeType l = str.length();
		if (l < 4) err("Pvector " + str + " string too short\n");
		SizeType last = l - 1;
		if (str_[last] != '>')
			weight_ = findWeightAndStripIt(str_);
	}

	void setString(PsimagLite::String newstring) { str_ = newstring; }

	void multiplyWeight(const RealType& factor) { weight_*= factor; }

	const PsimagLite::String& toString() { return str_; }

	const RealType& weight() const { return weight_; }

private:

	static RealType findWeightAndStripIt(PsimagLite::String str)
	{
		const SizeType l = str.length();
		if (l < 4) err("Pvector " + str + " string too short\n");
		PsimagLite::String buffer("");
		for (SizeType i = 0; i < l; ++i) {
			const SizeType j = l - i - 1;
			const unsigned char letter = str[j];
			if (letter == '*') break;
			if (!isAdigit(letter) && letter != '.' && letter != '+' && letter != '-')
				err("Wrong weight for vector " + str + "\n");

			buffer += letter;
		}

		return atoi(buffer.c_str());
	}

	static bool isAdigit(unsigned char letter)
	{
		return (letter > 47 && letter < 58);
	}

	PsimagLite::String str_;
	RealType weight_;
};
}
#endif // PVECTOR_H
