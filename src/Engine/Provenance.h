#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "Version.h" // do not commit this file, created dynamically

class Provenance {

public:


}; // Provenance

std::ostream& operator<<(std::ostream& os,const Provenance &prov)
{
	os<<"DMRG++: revision: "<<dmrgppRevision<<"\n";
	os<<"DMRG++: revision: "<<dmrgppDiff<<"\n";
	os<<"PsimagLite: revision: "<<psimagLiteRevision<<"\n";
	os<<"PsimagLite: diff: "<<psimagLiteDiff<<"\n";
	return os;
}

#endif // PROVENANCE_H
