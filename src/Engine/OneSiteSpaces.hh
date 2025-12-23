#ifndef ONESITESPACES_HH
#define ONESITESPACES_HH
#include "ProgramGlobals.h"

namespace Dmrg {

template <typename ModelType> class OneSiteSpaces {

public:

	OneSiteSpaces(SizeType site, ProgramGlobals::DirectionEnum dir, const ModelType& model)
	    : dir_(dir)
	{
		int siteAux = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM
		               || dir == ProgramGlobals::DirectionEnum::INFINITE)
		    ? site + 1
		    : site - 1;
		assert(siteAux >= 0);
		mainHilbert_ = model.hilbertSize(site);
		auxHilbert_ = model.hilbertSize(siteAux);
		zeroHilbert_ = model.hilbertSize(0);
	}

	void setDir(ProgramGlobals::DirectionEnum dir) { dir_ = dir; }

	SizeType hilbertMain() const { return mainHilbert_; }

	SizeType hilbertAux() const { return auxHilbert_; }

	SizeType hilbertZero() const { return zeroHilbert_; }

	ProgramGlobals::DirectionEnum direction() const { return dir_; }

private:

	ProgramGlobals::DirectionEnum dir_;
	SizeType mainHilbert_;
	SizeType auxHilbert_;
	SizeType zeroHilbert_;
};
}
#endif // ONESITESPACES_HH
