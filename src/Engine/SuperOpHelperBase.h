#ifndef SuperOpHelperBase_H
#define SuperOpHelperBase_H
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

class SuperOpHelperBase {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<bool, SizeType> PairBoolSizeType;

	SuperOpHelperBase(const VectorSizeType&,
	                  const VectorSizeType&,
	                  ProgramGlobals::DirectionEnum dir)
	    : dir_(dir) {}

	virtual ~SuperOpHelperBase() {}

	// This below is for a plaquette, and will have to be
	// written somewhere else
	// testing devel FIXME TODO
	virtual SizeType size() const { return 0; }

	virtual PairBoolSizeType leftOperatorIndex(SizeType) const
	{
		return PairBoolSizeType(false, 0);
	}

	virtual PairBoolSizeType rightOperatorIndex(SizeType) const
	{
		return PairBoolSizeType(false, 0);
	}

	const ProgramGlobals::DirectionEnum& dir() const { return dir_; }

private:

	ProgramGlobals::DirectionEnum dir_;
};
}
#endif // SuperOpHelperBase_H
