#ifndef SuperOpHelperPlaquette_H
#define SuperOpHelperPlaquette_H
#include "ProgramGlobals.h"
#include "SuperOpHelperBase.h"
#include "Vector.h"
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

template <typename SuperGeometryType, typename ParamsType>
class SuperOpHelperPlaquette : public SuperOpHelperBase<SuperGeometryType, ParamsType> {

public:

	typedef SuperOpHelperBase<SuperGeometryType, ParamsType> BaseType;
	typedef typename BaseType::VectorSizeType                VectorSizeType;
	typedef typename BaseType::PairBoolSizeType              PairBoolSizeType;
	using PairMetaOpForConnection = typename BaseType::PairMetaOpForConnection;

	SuperOpHelperPlaquette(const SuperGeometryType& superGeometry)
	    : BaseType(superGeometry)
	    , smaxOrEmin_(0)
	    , newSite_(0)
	{ }

	void setToProduct(SizeType smaxOrEmin, SizeType newSite, ProgramGlobals::DirectionEnum dir)
	{
		smaxOrEmin_ = smaxOrEmin;
		newSite_    = newSite;
		BaseType::setToProduct(smaxOrEmin, newSite, dir);
		isTriangle_ = (newSite & 1);
		if (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			isTriangle_ = !isTriangle_;
	}

	// This below is for a plaquette, and will have to be
	// written somewhere else
	// testing devel FIXME TODO
	SizeType size() const { return 1; }

	PairMetaOpForConnection finalIndices(const VectorSizeType&          hItems,
	                                     ProgramGlobals::ConnectionEnum type,
	                                     SizeType                       rightBlockSize) const
	{
		assert(hItems.size() == 4);
		constexpr int NON_LOCAL = -1;
		if (smaxOrEmin_ == 0) {
			// 0 x (1, 2, 3)
			MetaOpForConnection left { static_cast<int>(hItems[0]), 0, 'N' };
			MetaOpForConnection right { NON_LOCAL, encodeNonLocalEnv(1, 3), 'N' };
			return PairMetaOpForConnection(left, right);
		}

		SizeType nsites = BaseType::superGeometry().numberOfSites();

		if (smaxOrEmin_ + 2 == nsites) {
			// (n - 4, n - 3, n - 2) x (n - 1)
			MetaOpForConnection left { NON_LOCAL,
				                   encodeNonLocalSys(nsites - 4, 3),
				                   'N' };
			MetaOpForConnection right { static_cast<int>(hItems[3]), 0, 'N' };
			return PairMetaOpForConnection(left, right);
		}

		if (smaxOrEmin_ & 1) {
			// (s - 1, s) x (s + 1, s + 2)
			MetaOpForConnection left { NON_LOCAL,
				                   encodeNonLocalSys(smaxOrEmin_ - 1, 2),
				                   'N' };
			SizeType            start = rightBlockSize + smaxOrEmin_ - 2;
			MetaOpForConnection right { NON_LOCAL, encodeNonLocalEnv(start, 2), 'N' };
			return PairMetaOpForConnection(left, right);
		}

		// (s - 2, s - 1, s) x (s + 1)
		// return PairSizeType(encodeNonLocalSys(smaxOrEmin_ - 2, 3), smaxOrEmin_ + 1);

		// (s) x (s + 1, s + 2, s + 3)
		// return PairSizeType(smaxOrEmin_, encodeNonLocalEnv(smaxOrEmin_ + 1, 3));
		throw PsimagLite::RuntimeError("How to encode two plaquettes: depends on hItems\n");
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

	/*
	 *  (0, 1) --> 0
	 *  (0, 1, 2) --> 1
	 *  (2, 3) --> 2
	 *  (2, 3, 4) --> 3
	 *  (4, 5) --> 4
	 *  etc.
	 */
	SizeType encodeNonLocalSys(SizeType start, SizeType size) const
	{
		assert(size == 2 || size == 3);
		assert((start + 1) & 1);
		assert(start + size <= smaxOrEmin_ + 1);
		return (size == 2) ? start : start + 1;
	}

	/*
	 *  (n-2, n-1) --> 0
	 *  (n-3, n-2, n-1) --> 1
	 *  (n-4, n-3) --> 2
	 *  (n-5, n-4, n-3) --> 3
	 *  etc.
	 */
	SizeType encodeNonLocalEnv(SizeType start, SizeType size) const
	{
		SizeType nsites = BaseType::superGeometry().numberOfSites();
		assert(size == 2 || size == 3);
		assert(start > smaxOrEmin_);
		assert(start + size <= nsites);
		assert(nsites >= start + 2);
		return nsites - start - 2;
	}

	SizeType smaxOrEmin_;
	SizeType newSite_;
	bool     isTriangle_;
};
}
#endif // SuperOpHelperPlaquette_H
