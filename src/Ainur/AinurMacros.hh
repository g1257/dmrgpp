#ifndef AINURMACROS_HH
#define AINURMACROS_HH
#include <string>
#include "../Vector.h"

namespace PsimagLite {

class AinurMacros {

public:

	struct NativeMacro {
		std::string type;
		std::string name;
		std::string value;
	};

	AinurMacros()
	{
		nativeMacros_.push_back({"function", "AinurFromFile", "!AinurFromFile"});
	}

	SizeType total() const { return nativeMacros_.size(); }

	const NativeMacro& nativeMacro(SizeType ind) const
	{
		assert(ind < nativeMacros_.size());
		return nativeMacros_[ind];
	}

	std::string procNativeMacro(const std::string& line)
	{
		throw RuntimeError("procNativeMacro not implemented\n");
	}

private:

	std::vector<NativeMacro> nativeMacros_;
};

}
#endif // AINURMACROS_HH
