#ifndef PVECTOR_H
#define PVECTOR_H
#include "Vector.h"
#include <cstdlib>

namespace Dmrg {

template<typename ComplexOrRealType>
class Pvector {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

	// |P0>=(c?0[0]'*c?0[1]' +  c?1[0]'*c?1[1] - c?0[1]'*c?0[0] - c?1[1]'*c?1[0])|gs>*weight
	// The weight is optional
	Pvector(PsimagLite::String str)
	    : weight_(1.0)
	{
		// find the weight first
		SizeType l = str.length();
		if (l < 4) err("Pvector " + str + " string too short\n");
		SizeType last = l - 1;
		if (str[last] != '>')
			weight_ = findWeightAndStripIt(str);
		vStr_.push_back(str);
	}

	void setAsDone() {vStr_.push_back("DONE"); }

	void pushString(PsimagLite::String newstring) { vStr_.push_back(newstring); }

	void multiplyWeight(const RealType& factor) { weight_*= factor; }

	void sum(const Pvector& other, PsimagLite::String str)
	{
		const SizeType n = vStr_.size();
		if (n == 0 || vStr_[n - 1] != "DONE" || other.vStr_.size() == 0)
			err("Pvector::sum\n");

		PsimagLite::String def = vStr_[0] + other.vStr_[0];
		vStr_.clear();
		vStr_.resize(3);
		vStr_[0] = def;
		vStr_[1] = str;
		vStr_[2] = "DONE";
		// weights??
 	}

	const PsimagLite::String& firstName() const
	{
		if (vStr_.size() == 0)
			err("Pvector has no name\n");
		return vStr_[0];
	}

	const PsimagLite::String& lastName() const
	{
		const SizeType n = vStr_.size();
		if (n == 0)
			err("Pvector has no name\n");
		return vStr_[n - 1];
	}

	bool hasAnyName(PsimagLite::String str) const
	{
		const SizeType n = vStr_.size();
		for (SizeType i = 0; i < n; ++i)
			if (vStr_[i] == str) return true;
		return false;
	}

	const RealType& weight() const { return weight_; }

	SizeType size() const { return vStr_.size(); }

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

	VectorStringType vStr_;
	RealType weight_;
};
}
#endif // PVECTOR_H
