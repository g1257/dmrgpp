#ifndef GETBRAORKET_H
#define GETBRAORKET_H
#include "Vector.h"
#include <cstdlib>

namespace PsimagLite {

class GetBraOrKet {

public:

	GetBraOrKet(String braOrKet) : braOrKet_(braOrKet)
	{
		const SizeType l = braOrKet.length();
		if (l < 2)
			err("GetBraOrKet " + braOrKet + "\n");
		const SizeType last = l - 1;
		if (braOrKet_[0] == '|' && braOrKet_[last] == '>' && l > 2)
			braOrKet_ = braOrKet.substr(1, l - 2);
		if (braOrKet_[0] == '<' && braOrKet_[last] == '|' && l > 2)
			braOrKet_ = braOrKet.substr(1, l - 2);
	}

	SizeType operator()() const
	{
		if (braOrKet_ == "gs")
			return 0;

		int ind = getPtype(braOrKet_);
		if (ind <= 0)
			err("Malformed braket " + braOrKet_ + "\n");
		return ind;
	}

	static int getPtype(String str)
	{
		// str == P\d+
		if (str.length() < 2) return -1;
		if (str[0] != 'P') return -1;
		String number("");
		for (SizeType i = 1; i < str.length(); ++i) {
			number += str[i];
			unsigned char x = str[i];
			if (x < 48 || x > 57) return -1;
		}

		return atoi(number.c_str()) + 1;
	}

private:

	String braOrKet_;
};
}
#endif // GETBRAORKET_H
