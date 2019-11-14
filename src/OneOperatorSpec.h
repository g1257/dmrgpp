#ifndef ONEOPERATORSPEC_H
#define ONEOPERATORSPEC_H
#include "Vector.h"
#include <cstdlib>

namespace PsimagLite {
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

	struct SiteSplit {

		SiteSplit(bool hasSiteString_, String root_, String siteString_)
		    : hasSiteString(hasSiteString_), root(root_), siteString(siteString_)
		{}

		bool hasSiteString;
		String root;
		String siteString;
	};

	static SiteSplit extractSiteIfAny(PsimagLite::String name)
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

		if (firstIndex < 0 && lastIndex < 0) return SiteSplit(false, name, "");

		bool b1 = (firstIndex < 0 && lastIndex >= 0);
		bool b2 = (firstIndex >= 0 && lastIndex < 0);
		if (b1 || b2)
			err(name + " has unmatched [ or ]\n");

		String str = name.substr(0, firstIndex);
		str += name.substr(lastIndex + 1, name.length() - lastIndex);
		String siteString = name.substr(firstIndex+1,lastIndex-1).c_str();
		return SiteSplit(true, str, siteString);
	}

	static bool isNonNegativeInteger(const String& s)
	{
		return !s.empty() &&
		        std::find_if(s.begin(),
		                     s.end(),
		                     [](char c) { return !std::isdigit(c); }) == s.end();
	}

	static SizeType strToNumberOfFail(String s)
	{
		if (!isNonNegativeInteger(s))
			err("string " + s + " is not a NonNegativeInteger\n");
		return atoi(s.c_str());
	}

	SizeType dof;
	String label;
	bool transpose;
}; // struct OneOperatorSpec

}
#endif // ONEOPERATORSPEC_H
