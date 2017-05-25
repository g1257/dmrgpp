#ifndef PSI_PSIMAGLITE_H
#define PSI_PSIMAGLITE_H

#include <iostream>
#include <utility>
#include "Concurrency.h"
#include "AnsiColors.h"

namespace PsimagLite {

std::ostream& operator<<(std::ostream&,const std::pair<SizeType,SizeType>&);

std::istream& operator>>(std::istream&,std::pair<SizeType,SizeType>&);

SizeType log2Integer(SizeType x);

void err(String);

class PsiApp {
public:

	PsiApp(String appName, int* argc, char*** argv, int nthreads)
	    : concurrency_(argc,argv,nthreads)
	{
		chekSizeType();
	}

private:

	void chekSizeType()
	{
		if (sizeof(SizeType) == libSizeOfSizeType_) return;
		std::string msg("PsimagLite compiled with -DUSE_LONG but");
		msg += "application without. Or viceversa.\n";
		throw std::runtime_error(msg);
	}

	static const int libSizeOfSizeType_;

	Concurrency concurrency_;
};
} // namespace PsimagLite

void err(PsimagLite::String);

#endif // PSI_PSIMAGLITE_H

