#ifndef OPTIONS_WRITEABLE_H
#define OPTIONS_WRITEABLE_H

#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include "Split.h"

namespace PsimagLite {

class OptionsWriteable {

public:

	class OptionsReadable {
	public:
		OptionsReadable(OptionsWriteable& optsWrite,const std::string& optsString)
		{
			optsWrite.set(optsThatAreSet_,optsString);
		}

		bool isSet(const std::string& thisOption) const
		{
			bool b = (find(optsThatAreSet_.begin(),optsThatAreSet_.end(),thisOption)==optsThatAreSet_.end());
			return (!b);
		}
	private:
		std::vector<std::string> optsThatAreSet_;

	}; // class OptionsReadable

	OptionsWriteable(std::vector<std::string>& registeredOptions)
		: registeredOptions_(registeredOptions)
	{}

	void set(std::vector<std::string>& optsThatAreSet,const std::string& opts)
	{
		split(optsThatAreSet,opts.c_str(),',');
		for (size_t i=0;i<optsThatAreSet.size();i++) {
			bool b = (find(registeredOptions_.begin(),registeredOptions_.end(),optsThatAreSet[i])==registeredOptions_.end());
			if (b) {
				std::string s(__FILE__);
				s += ": Unknown option " + optsThatAreSet[i] + "\n";
				throw std::runtime_error(s.c_str());
			}
		}
	}

private:
	std::vector<std::string> registeredOptions_;

}; // class OptionsWriteable

} // namespace PsimagLite
#endif // OPTIONS_WRITEABLE_H
