#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "../../PsimagLite/src/Version.h"
#include "../Version.h"
#include "AllocatorCpu.h"
#include <iostream>
#include "AnsiColors.h"
#include <sstream>
#include "../GitRevision.h"

class Provenance {

public:

	static PsimagLite::String compiledMicroArch()
	{
#ifndef MICRO_ARCH
#error "Please run ./configure.pl ...\n";
		return "";
#else
		return MICRO_ARCH;
#endif
	}

	static PsimagLite::String logo(PsimagLite::String appName)
	{
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<appName<<"\x1b[38;5;240m";
		msg<<" [master "<<DMRGPP_VERSION<<"] "<<PsimagLite::AnsiColor::reset;
		PsimagLite::String ctOpts("");
#ifdef USE_LONG
		ctOpts += "LONG ";
#endif
#ifdef USE_FLOAT
		ctOpts += "FLOAT ";
#endif
#ifdef USE_SIGNALS
		ctOpts += "SIGNALS ";
#endif
#ifdef USE_GSL
		ctOpts += "GSL ";
#endif
#ifdef USE_BOOST
		ctOpts += "BOOST ";
#endif
#ifdef PLUGIN_SC
		ctOpts += "PLUGIN_SC ";
#endif
#ifndef NDEBUG
		ctOpts += "DEBUG ";
#endif

		if (ctOpts != "")
			msg<<"["<<ctOpts<<"]";

		return msg.str();
	}
}; // Provenance

std::ostream& operator<<(std::ostream& os,const Provenance&);

#endif // PROVENANCE_H

