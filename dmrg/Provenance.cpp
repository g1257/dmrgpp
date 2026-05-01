#include "Provenance.h"

std::ostream& operator<<(std::ostream& os, const Provenance&)
{
	os << "DMRG++ version " << DMRGPP_VERSION << " " << DMRGPP_GIT_REV << "\n";
	return os;
}
