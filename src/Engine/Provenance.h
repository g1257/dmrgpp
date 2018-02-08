#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "../../PsimagLite/src/Version.h"
#include "../Version.h"
#include "AllocatorCpu.h"
#include <iostream>
#include "AnsiColors.h"
#include <sstream>

class Provenance {

public:

	static PsimagLite::String logo(PsimagLite::String appName)
	{
		PsimagLite::OstringStream msg;
		msg<<appName<<"\x1b[38;5;240m";
		msg<<" [master "<<DMRGPP_VERSION<<"] "<<PsimagLite::AnsiColor::reset;
#ifdef PLUGIN_SC
		msg<<"[PLUGIN_SC] ";
#endif

		return msg.str();
	}
}; // Provenance

std::ostream& operator<<(std::ostream& os,const Provenance&);

#endif // PROVENANCE_H

