#ifndef PSI_PSIMAGLITE_H
#define PSI_PSIMAGLITE_H

#include <iostream>
#include <utility>
#include "Vector.h"

namespace PsimagLite {

std::ostream& operator<<(std::ostream&,const std::pair<SizeType,SizeType>&);

std::istream& operator>>(std::istream&,std::pair<SizeType,SizeType>&);

class PsiApp {
public:

	PsiApp(String appName)
	{
		chekSizeType();
	}

private:

	void chekSizeType()
	{
		if (sizeof(SizeType) == libSizeOfSizeType()) return;
		std::string msg("PsimagLite compiled with -DUSE_LONG but");
		msg += "application without. Or viceversa.\n";
		throw std::runtime_error(msg);
	}

	int libSizeOfSizeType();
};

} // namespace PsimagLite

#endif // PSI_PSIMAGLITE_H

