#ifndef SuperOpHelperBase_H
#define SuperOpHelperBase_H
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

template<typename SuperGeometryType, typename ParametersType>
class SuperOpHelperBase {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<bool, SizeType> PairBoolSizeType;

	SuperOpHelperBase(const SuperGeometryType& superGeometry)
	    : superGeometry_(superGeometry)
	{}

	virtual ~SuperOpHelperBase() {}

	virtual void setToProduct(SizeType, SizeType, ProgramGlobals::DirectionEnum dir)
	{
		dir_ = dir;
	}

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

	virtual SizeType leftIndex(VectorSizeType&, SizeType) const
	{
		throw PsimagLite::RuntimeError("SuperOpHelperBase::leftIndex\n");
	}

	virtual SizeType rightIndex(VectorSizeType&, SizeType) const
	{
		throw PsimagLite::RuntimeError("SuperOpHelperBase::rightIndex\n");
	}

	// non virtual below

	const SuperGeometryType& superGeometry() const
	{
		return superGeometry_;
	}

	const ProgramGlobals::DirectionEnum& dir() const { return dir_; }

private:

	const SuperGeometryType& superGeometry_;
	ProgramGlobals::DirectionEnum dir_;
};
}
#endif // SuperOpHelperBase_H
