#ifndef SuperOpHelperBase_H
#define SuperOpHelperBase_H
#include "MetaOpForConnection.hh"
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

template <typename SuperGeometryType, typename ParametersType> class SuperOpHelperBase {

public:

	using VectorSizeType          = PsimagLite::Vector<SizeType>::Type;
	using PairBoolSizeType        = std::pair<bool, SizeType>;
	using PairMetaOpForConnection = std::pair<MetaOpForConnection, MetaOpForConnection>;

	SuperOpHelperBase(const SuperGeometryType& superGeometry)
	    : superGeometry_(superGeometry)
	{ }

	virtual ~SuperOpHelperBase() { }

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

	virtual PairMetaOpForConnection
	finalIndices(const VectorSizeType&, ProgramGlobals::ConnectionEnum, SizeType) const
	{
		throw PsimagLite::RuntimeError("SuperOpHelperBase::finalIndices4sites\n");
	}

	// non virtual below

	const SuperGeometryType& superGeometry() const { return superGeometry_; }

	const ProgramGlobals::DirectionEnum& dir() const { return dir_; }

private:

	const SuperGeometryType&      superGeometry_;
	ProgramGlobals::DirectionEnum dir_;
};
}
#endif // SuperOpHelperBase_H
