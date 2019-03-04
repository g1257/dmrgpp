#ifndef ONEOPERATORSPEC_H
#define ONEOPERATORSPEC_H
#include "Vector.h"
#include <cstdlib>

namespace Dmrg {
struct OneOperatorSpec {

	OneOperatorSpec(PsimagLite::String label_)
	    : dof(0),label(label_),transpose(false)
	{
		SizeType lastIndex = label.length();
		if (lastIndex > 0) lastIndex--;
		if (label[lastIndex] == '\'') {
			label = label.substr(0,lastIndex);
			transpose = true;
		}

		label_ = label;

		SizeType i = 0;
		for (; i < label.length(); ++i) {
			if (label[i] == '?') break;
		}

		if (i == label.length()) return;

		if (i + 1 == label.length())
			err("WRONG op. spec. " + label_ + ", nothing after ?\n");

		label = label_.substr(0, i);
		dof = atoi(label_.substr(i + 1, label_.length()).c_str());
	}

	static int extractSiteIfAny(PsimagLite::String& name)
	{
		int firstIndex = -1;
		int lastIndex = -1;
		for (SizeType i = 0; i < name.length(); ++i) {
			if (name[i] == '[') {
				firstIndex = i;
				continue;
			}

			if (name[i] == ']') {
				lastIndex = i;
				continue;
			}
		}

		if (firstIndex < 0 && lastIndex < 0) return -1;

		bool b1 = (firstIndex < 0 && lastIndex >= 0);
		bool b2 = (firstIndex >= 0 && lastIndex < 0);
		if (b1 || b2) {
			PsimagLite::String str("Braket operator ");
			err(name + " has unmatched [ or ]\n");
		}

		PsimagLite::String str = name.substr(0, firstIndex);
		str += name.substr(lastIndex + 1, name.length() - lastIndex);
		int site = atoi(name.substr(firstIndex+1,lastIndex-1).c_str());
		name = str;
		return site;
	}

	SizeType dof;
	PsimagLite::String label;
	bool transpose;
}; // struct OneOperatorSpec

}
#endif // ONEOPERATORSPEC_H
