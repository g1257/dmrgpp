#ifndef METAOPFORCONNECTION_HH
#define METAOPFORCONNECTION_HH
#include "AllocatorCpu.h"

namespace Dmrg {

struct MetaOpForConnection {

	MetaOpForConnection(int site1, SizeType index1, char modifier1)
	    : site(site1), index(index1), modifier(modifier1)
	{}

	int site; // -1 means non-local
	SizeType index;
	char modifier;
};

}
#endif // METAOPFORCONNECTION_HH
