#ifndef RUNFINISHED_H
#define RUNFINISHED_H
#include "PsimagLite.h"

namespace Dmrg {

class RunFinished {

public:

	RunFinished(bool enabled)
	    : enabled_(enabled)
	    , checked_(false)
	    , value_(false)
	{
	}

	bool OK(PsimagLite::String filename)
	{
		if (!enabled_)
			return false;
		if (!checked_)
			checkRun(filename);
		return value_;
	}

	void printTermination(PsimagLite::String filename)
	{
		static const PsimagLite::String str = "File " + filename + " exists, " + "and you chose no clobber. Refusing to run\n";
		std::cerr << str;
		std::cout << str;
	}

private:

	void checkRun(PsimagLite::String filename)
	{
		checked_ = true;
		value_ = false;
		std::ifstream fin(filename.c_str(), std::ios::ate);
		if (!fin || fin.bad() || !fin.good()) {
			value_ = false;
			return;
		}

		char schar1[1024];
		std::streampos size = fin.tellg();
		int i = 1;
		while (i <= size) {
			fin.seekg(-i, std::ios::end);
			fin.getline(schar1, 1023);
			PsimagLite::String str(schar1);
			if (str.find("Turning off the engine.") != PsimagLite::String::npos) {
				value_ = true;
				break;
			}
			i++;
		}
		fin.close();
	}

	bool enabled_;
	bool checked_;
	bool value_;
};
}
#endif // RUNFINISHED_H
