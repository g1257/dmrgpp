#ifndef GETBRAORKET_H
#define GETBRAORKET_H
#include "Vector.h"
#include <cstdlib>

namespace PsimagLite {

class GetBraOrKet {

public:

	enum class Kind {E, P};

	typedef std::pair<SizeType, SizeType> PairSizeType;

	GetBraOrKet(String braOrKet)
	    : isKet_(true), braOrKet_(braOrKet), kind_(Kind::E), pair_(0, 0)
	{
		const SizeType l = braOrKet.length();
		if (l < 2)
			err("GetBraOrKet " + braOrKet + "\n");
		const SizeType last = l - 1;
		if (braOrKet_[0] == '|' && braOrKet_[last] == '>' && l > 2)
			braOrKet_ = braOrKet.substr(1, l - 2);
		if (braOrKet_[0] == '<' && braOrKet_[last] == '|' && l > 2) {
			braOrKet_ = braOrKet.substr(1, l - 2);
			isKet_ = false;
		}

		getKind(braOrKet_);
	}

	bool isKet() const { return isKet_; }

	const String& toString() const { return braOrKet_; }

	bool isPvector() const { return (kind_ == Kind::P); }

	SizeType levelIndex() const
	{
		checkIfPvector(false);
		return pair_.second;
	}

	SizeType sectorIndex() const
	{
		checkIfPvector(false);
		return pair_.first;
	}

	SizeType pIndex() const
	{
		checkIfPvector(true);
		return pair_.first;
	}

private:

	void checkIfPvector(bool b) const
	{
		bool isP = (kind_ == Kind::P);
		if (isP == b) return;
		throw PsimagLite::RuntimeError("Internal ERROR: checkIfPpvector\n");
	}

	void getKind(String str)
	{
		if (str.length() < 2)
			err("GetBraOrKet:: " + str + "too short\n");

		if (str[0] == 'P') {
			kind_ = Kind::P;
			pair_.first = getNumberFrom(str, 1); // modifies str
		} else if (str == "gs") {
			kind_ = Kind::E;
			return;
		} else if (str == "time") { // legacy name
			kind_ = Kind::P;
			return;
		} else if (str[0] == 'Q') {
			kind_ = Kind::E;
			pair_.first = getNumberFrom(str, 1); // modifies str
		}

		if (str.size() < 2)
			return;

		if (str[0] == 'X') {
			pair_.second = getNumberFrom(str, 1); // modifies str
		} else {
			err("A vector spec can only start with P, gs, X or Q " + str + "\n");
		}
	}

	static SizeType getNumberFrom(String& str, SizeType start)
	{
		String number("");
		SizeType i = start;
		for (; i < str.length(); ++i) {
			unsigned char x = str[i];
			if (x < 48 || x > 57) break;
			number += str[i];
		}

		if (number == "")
			RuntimeError("getNumberFrom: no number after character location " + ttos(start) +
			             " for string " + str + "\n");

		str = str.substr(i, str.length() - i);

		return atoi(number.c_str());
	}

	bool isKet_;
	String braOrKet_;
	Kind kind_;
	PairSizeType pair_;
};
}
#endif // GETBRAORKET_H
