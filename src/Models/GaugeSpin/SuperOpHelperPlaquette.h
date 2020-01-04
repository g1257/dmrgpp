#ifndef SuperOpHelperPlaquette_H
#define SuperOpHelperPlaquette_H
#include "ProgramGlobals.h"
#include "Vector.h"
#include "SuperOpHelperBase.h"

namespace Dmrg {

template<typename SuperGeometryType, typename ParamsType>
class SuperOpHelperPlaquette : public SuperOpHelperBase<SuperGeometryType, ParamsType> {

public:

	typedef SuperOpHelperBase<SuperGeometryType, ParamsType> BaseType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::PairBoolSizeType PairBoolSizeType;

	SuperOpHelperPlaquette(const SuperGeometryType& superGeometry)
	    : BaseType(superGeometry), smaxOrEmin_(0), newSite_(0)
	{}

	void setToProduct(SizeType smaxOrEmin,
	                  SizeType newSite,
	                  ProgramGlobals::DirectionEnum dir)
	{
		smaxOrEmin_ = smaxOrEmin;
		newSite_ = newSite;
		BaseType::setToProduct(smaxOrEmin, newSite, dir);
		isTriangle_ = (newSite & 1);
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

	SizeType leftIndex(VectorSizeType& sysSites, SizeType) const
	{
		SizeType last = sysSites.size();
		assert(last > 0);
		--last;

	}

	SizeType rightIndex(VectorSizeType&, SizeType) const
	{
		throw PsimagLite::RuntimeError("SuperOpHelperBase::rightIndex\n");
	}

private:

	SizeType smaxOrEmin_;
	SizeType newSite_;
	bool isTriangle_;
};
}
#endif // SuperOpHelperPlaquette_H
