#ifndef SuperOpHelperPlaquette_H
#define SuperOpHelperPlaquette_H
#include "ProgramGlobals.h"
#include "Vector.h"
#include "SuperOpHelperBase.h"
/*
 * Non local ops.
 *
 * Only even sites have non-local ops, and only one each.
 *
 * 0 --> 2 --> 4 --> ...
 * |     |     |
 *
 * 1 --> 3 --> 5 --> ...
 *
 * even sites have a non-local op. equal to the left top angle
 *
 * x --> x + 2
 * |
 *
 * x + 1
 *
 * These operator is called A.
 *
 * Therefore the plaquettes that need to be added as
 * connections between system and environment (for
 * a partition where smax is the maximum site of the
 * system and H_R is the addition to the right Hamiltonian)
 * are as follows.
 *
 * smax == 0 A_0 * sx_{3, y} * sx_{1, x}, H_R = 0
 *
 * smax == 1 same as smax==0, H_R = 0
 *
 * smax == 2 has two contributions
 *  same as smax==0 plus
 *  A_2 * sx_{5, y} * sx_{3, x}, H_R = 0
 *
 * smax == 3 here we add plaquette 0-2-3-1 to H_R
 * and the tensor product to be added is
 * A_2 * sx_{5, y} * sx_{3, x}
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

	PairSizeType finalIndices4sites(const VectorSizeType& hItems,
	                                ProgramGlobals::ConnectionEnum type) const
	{
		if (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) {
			return finalIndices4sitesSysEnv(hItems);
		} else if (type == ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM) {
			return finalIndices4sitesEnvSys(hItems);
		}

		throw PsimagLite::RuntimeError("Internal error, unexpected type\n");
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


	PairSizeType finalIndices4sitesSysEnv(const VectorSizeType& hItems) const
	{
		throw PsimagLite::RuntimeError("finalIndices4sitesSysEnv\n");
	}

	PairSizeType finalIndices4sitesEnvSys(const VectorSizeType& hItems) const
	{
		throw PsimagLite::RuntimeError("finalIndices4sitesEnvSys\n");
	}

	SizeType smaxOrEmin_;
	SizeType newSite_;
	bool isTriangle_;
};
}
#endif // SuperOpHelperPlaquette_H
