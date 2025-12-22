#ifndef LOOPSITEDIRECTION_HH
#define LOOPSITEDIRECTION_HH
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg
{

struct LoopSiteDirection {
	SizeType loopIndex;
	SizeType site;
	ProgramGlobals::DirectionEnum direction;
};

}
#endif // LOOPSITEDIRECTION_HH
