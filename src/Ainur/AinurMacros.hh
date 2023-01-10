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

	AinurMacros() : AINUR_FROM_FILE("AinurFromFile")
	{
		nativeMacros_.push_back({"function", AINUR_FROM_FILE, "!" + AINUR_FROM_FILE});
	}

	SizeType total() const { return nativeMacros_.size(); }

	const NativeMacro& nativeMacro(SizeType ind) const
	{
		assert(ind < nativeMacros_.size());
		return nativeMacros_[ind];
	}

	std::string procNativeMacro(const std::string& line)
	{
		if (line.substr(1, AINUR_FROM_FILE.size()) == AINUR_FROM_FILE) {
			return procAinurFromFile(line);
		}

		throw RuntimeError("Unknown native macro for " + line + "\n");
	}

private:

	std::string procAinurFromFile(const std::string& line)
	{
		SizeType start = AINUR_FROM_FILE.size() + 1;
		SizeType len = line.size();
		assert(len > start);
		SizeType rest = len - start;
		std::string content = line.substr(start, rest);
		throw RuntimeError("unimplemented for " + content + "\n");
	}

	const std::string AINUR_FROM_FILE;
	std::vector<NativeMacro> nativeMacros_;
};

}
#endif // AINURMACROS_HH
