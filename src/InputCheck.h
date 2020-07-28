/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef DMFT_INPUT_CHECK_H
#define DMFT_INPUT_CHECK_H
#include <vector>
#include <stdexcept>
#include "../../PsimagLite/src/Options.h"
#include "Geometry/Geometry.h"
//#include "ProgramGlobals.h"

namespace Dmft {

class InputCheck {

	typedef PsimagLite::Options::Readable OptionsReadableType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;

public:

	InputCheck() : optsReadable_(0)
	{
		//knownLabels_.push_back("TotalNumberOfSites");
	}

	~InputCheck()
	{
		if (optsReadable_!=0) delete optsReadable_;
	}

	PsimagLite::String import() const
	{
		//PsimagLite::String str = PsimagLite::Geometry<int,int,ProgramGlobals>::import();
		PsimagLite::String str;

		str += "integer FicticiousBeta;\n";
		str += "real ChemicalPotential;\n";
		str += "integer Matsubaras;\n";
		str += "integer NumberOfKpoints;\n";
		str += "integer NumberOfBathPoints;\n";
		str += "integer DmftNumberOfIterations;\n";
		str += "real DmftTolerance;\n";
		str += "real MinParamsDelta;\n";
		str += "real MinParamsDelta2;\n";
		str += "real MinParamsTolerance;\n";
		str += "integer MinParamsMaxIter;\n";
		str += "integer MinParamsVerbose;\n";

		return str;
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		return false;
	}

	bool check(const PsimagLite::String& label,
	           const PsimagLite::String& vec,
	           SizeType line) const
	{
		return false;
	}

	bool checkSimpleLabel(const PsimagLite::String& label,
	                      SizeType line) const
	{
		for (SizeType i = 0; i < knownLabels_.size(); ++i)
			if (knownLabels_[i] == label) return true;
		PsimagLite::String msg("WARNING: Unknown label " + label +"\n");
		std::cout<<msg;
		std::cerr<<msg;
		return false;
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

private:

	bool checkForVector(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() == 0) return false;
		SizeType n = atoi(vec[0].c_str());
		return (vec.size() == n+1);
	}

	bool checkForMatrix(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() < 2) return false;
		SizeType row = atoi(vec[0].c_str());
		SizeType col = atoi(vec[1].c_str());
		SizeType n = row*col;
		return (vec.size() == n+2);
	}

	bool error1(const PsimagLite::String& message,SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());

	}

	OptionsReadableType* optsReadable_;
	VectorStringType allowedFileOptions_;
	VectorStringType knownLabels_;
}; // class InputCheck
} // namespace Dmft

/*@}*/
#endif

