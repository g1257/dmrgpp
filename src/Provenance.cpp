#include "Provenance.h"
#include "Version.h"

std::ostream& operator<<(std::ostream& os,const Provenance&)
{
	os<<"DMRG++ version "<<DMRGPP_VERSION<<"\n";
	os<<"PsimagLite version "<<PSIMAGLITE_VERSION<<"\n";
	return os;
}

