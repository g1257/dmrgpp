#ifndef METAOPFORCONNECTION_HH
#define METAOPFORCONNECTION_HH
#include "AllocatorCpu.h"

namespace Dmrg {

struct MetaOpForConnection {
	int      site; // -1 means non-local
	SizeType index;
	char     modifier;
};

}
#endif // METAOPFORCONNECTION_HH
