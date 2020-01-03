#ifndef SuperOpHelperPlaquette_H
#define SuperOpHelperPlaquette_H
#include "ProgramGlobals.h"
#include "Vector.h"
#include "SuperOpHelperBase.h"

namespace Dmrg {

class SuperOpHelperPlaquette : public SuperOpHelperBase {

public:

	typedef SuperOpHelperBase BaseType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::PairBoolSizeType PairBoolSizeType;

	SuperOpHelperPlaquette(const VectorSizeType& bigBlock,
	                       const VectorSizeType& smallBlock,
	                       ProgramGlobals::DirectionEnum dir)
	    : BaseType(bigBlock, smallBlock, dir), bigBlock_(bigBlock), smallBlock_(smallBlock)
	{
		assert(smallBlock_.size() == 1);
		isTriangle_ = (smallBlock_[0] & 1);
		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			isTriangle_ = !isTriangle_;
	}

	// This below is for a plaquette, and will have to be
	// written somewhere else
	// testing devel FIXME TODO
	SizeType size() const { return 1; }

	PairBoolSizeType leftOperatorIndex(SizeType) const
	{
		if (isTriangle_) {
			if (BaseType::dir() == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				return PairBoolSizeType(true, 0);
			} else {
				return PairBoolSizeType(false, 0);
			}
		}

		return PairBoolSizeType(false, 0);
	}

	PairBoolSizeType rightOperatorIndex(SizeType) const
	{
		if (isTriangle_) {
			if (BaseType::dir() == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				return PairBoolSizeType(false, 0);
			} else {
				return PairBoolSizeType(false, 0);
			}
		}

		return PairBoolSizeType(false, 0);
	}

private:

	const VectorSizeType& bigBlock_;
	const VectorSizeType& smallBlock_;
	bool isTriangle_;
};
}
#endif // SuperOpHelperPlaquette_H
