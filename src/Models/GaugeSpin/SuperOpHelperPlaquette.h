#ifndef SuperOpHelperPlaquette_H
#define SuperOpHelperPlaquette_H
#include "ProgramGlobals.h"
#include "Vector.h"
#include "SuperOpHelperBase.h"
/*
 * Non local ops.
 *
 *
 * 0 --> 2 --> 4 --> ...
 * |     |     |
 *
 * 1 --> 3 --> 5 --> ...
 *
 * (0, 1)
 * (0, 1, 2)
 * (2, 3)
 * (2, 3, 4)
 * (4, 5)
 *
 *
 * etc.
 *
 */
namespace Dmrg {

template<typename SuperGeometryType, typename ParamsType>
class SuperOpHelperPlaquette : public SuperOpHelperBase<SuperGeometryType, ParamsType> {

public:

	typedef SuperOpHelperBase<SuperGeometryType, ParamsType> BaseType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename BaseType::PairBoolSizeType PairBoolSizeType;
	typedef typename BaseType::PairSizeType PairSizeType;

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

	PairSizeType finalIndices(const VectorSizeType& hItems,
	                          ProgramGlobals::ConnectionEnum type) const
	{
		if (smaxOrEmin_ == 0) {
			// 0 x (1, 2, 3)
			return PairSizeType(0, encodeNonLocal(1, 3));
		}

		SizeType nsites = BaseType::superGeometry().numberOfSites();

		if (smaxOrEmin_ + 2 == nsites) {
			// (n - 4, n - 3, n - 2) x (n - 1)
			return PairSizeType(encodeNonLocal(nsites - 4, 3), nsites - 1);
		}

		if (smaxOrEmin_ & 1) {
			// (s - 1, s) x (s + 1, s + 2)
			return PairSizeType(encodeNonLocal(smaxOrEmin_ - 1, 2),
			                    encodeNonLocal(smaxOrEmin_ + 1, 2));
		}

		// (s - 2, s - 1, s) x (s + 1)
		//return PairSizeType(encodeNonLocal(smaxOrEmin_ - 2, 3), smaxOrEmin_ + 1);

		// (s) x (s + 1, s + 2, s + 3)
		throw PsimagLite::RuntimeError("How to encode two plaquettes?\n");
	}

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

	static SizeType encodeNonLocal(SizeType start, SizeType size)
	{
		throw PsimagLite::RuntimeError("encodeNonLocal\n");
	}

	SizeType smaxOrEmin_;
	SizeType newSite_;
	bool isTriangle_;
};
}
#endif // SuperOpHelperPlaquette_H
