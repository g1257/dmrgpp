#include "Provenance.h"

namespace Dmft {
std::ostream& operator<<(std::ostream& os, const Provenance&)
{
	os<<"Cincuenta version "<<DMFT_VERSION<<" "<<DMFT_GIT_REV<<"\n";
	os<<"PsimagLite version "<<PSIMAGLITE_VERSION<<" "<<PSIMAGLITE_GIT_REV<<"\n";
	return os;
}
}

