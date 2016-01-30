#include "Provenance.h"
#include "Version.h" // do not commit this file, created dynamically

std::ostream& operator<<(std::ostream& os,const Provenance&)
{
	os<<"DMRG++: revision: "<<dmrgppRevision<<"\n";
	os<<"DMRG++: diff: "<<dmrgppDiff<<"\n";
	os<<"PsimagLite: revision: "<<psimagLiteRevision<<"\n";
	os<<"PsimagLite: diff: "<<psimagLiteDiff<<"\n";
	return os;
}

