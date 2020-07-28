#include "Provenance.h"

std::ostream& operator<<(std::ostream& os, const Provenance&)
{
	os<<"Dmft version "<<DMFT_VERSION<<" "<<DMFT_GIT_REV<<"\n";
	os<<"PsimagLite version "<<PSIMAGLITE_VERSION<<" "<<PSIMAGLITE_GIT_REV<<"\n";
	return os;
}

