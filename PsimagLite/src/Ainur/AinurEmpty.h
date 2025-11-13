#ifndef _AINUR_EMPTY_H_
#define _AINUR_EMPTY_H_
#include "PsimagLite.h"
#include "TypeToString.h"
#include "Vector.h"
#include <fstream>
#include <iostream>

namespace PsimagLite {

class Ainur {

public:

	Ainur(String str)
	    : dummy_("")
	{
		errorMessage();
	}

	String& prefix() { return dummy_; }

	const String& prefix() const { return dummy_; }

	void printUnused(std::ostream& os) const { errorMessage(); }

	void printAll(std::ostream& os) const { errorMessage(); }

	template <typename SomeType>
	void readValue(SomeType& t, String label) const
	{
		errorMessage();
	}

	std::string resolve(const std::string& str) const
	{
		errorMessage();
		return "";
	}

	template <typename SomeMapType>
	void setMap(SomeMapType&) const
	{
		errorMessage();
	}

private:

	void errorMessage() const
	{
		err("To use Ainur, you need boost-devel, and compile with "
		    "-DUSE_BOOST\n");
	}

	String dummy_;
};

} // namespace PsimagLite
#endif // _AINUR_EMPTY_H_
