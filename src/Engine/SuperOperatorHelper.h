#ifndef SUPEROPERATORHELPER_H
#define SUPEROPERATORHELPER_H
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

template<typename ModelType>
class SuperOperatorHelper {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<bool, SizeType> PairBoolSizeType;

	SuperOperatorHelper(const ModelType& model,
	                    const VectorSizeType& bigBlock,
	                    const VectorSizeType& smallBlock,
	                    ProgramGlobals::DirectionEnum dir)
	    : model_(model), bigBlock_(bigBlock), smallBlock_(smallBlock), dir_(dir)
	{
		assert(smallBlock_.size() == 1);
		isTriangle_ = (smallBlock_[0] & 1);
		if (dir_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			isTriangle_ = !isTriangle_;
	}

	const ProgramGlobals::DirectionEnum& dir() const { return dir_; }

	// This below is for a plaquette, and will have to be
	// written somewhere else
	// testing devel FIXME TODO
	SizeType size() const { return 1; }

	PairBoolSizeType leftOperatorIndex(SizeType) const
	{
		if (isTriangle_) {
			if (dir_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
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
			if (dir_ == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				return PairBoolSizeType(false, 0);
			} else {
				return PairBoolSizeType(false, 0);
			}
		}

		return PairBoolSizeType(false, 0);
	}

private:

	const ModelType& model_;
	const VectorSizeType& bigBlock_;
	const VectorSizeType& smallBlock_;
	ProgramGlobals::DirectionEnum dir_;
	bool isTriangle_;
};
}
#endif // SUPEROPERATORHELPER_H
