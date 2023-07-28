#ifndef PROVENANCE_H
#define PROVENANCE_H
#include "../../PsimagLite/src/Version.h"
#include "../GitRevision.h"
#include "../Version.h"
#include "AllocatorCpu.h"
#include "AnsiColors.h"
#include "MatrixVectorKron/BatchedGemmInclude.hh"
#include <iostream>
#include <sstream>

class Provenance
{

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
		msg << appName << "\x1b[38;5;120";
		msg << " [features " << DMRGPP_VERSION << "] " << PsimagLite::AnsiColor::reset;
		PsimagLite::String ctOpts(Dmrg::BatchedGemmInclude::info());
#ifdef USE_SHORT
		ctOpts += " SHORT ";
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
#ifndef NDEBUG
		ctOpts += "DEBUG ";
#endif

		if (ctOpts != "")
			msg << "[" << ctOpts << "]";

		return msg.str();
	}
}; // Provenance

std::ostream& operator<<(std::ostream& os, const Provenance&);

#endif // PROVENANCE_H
