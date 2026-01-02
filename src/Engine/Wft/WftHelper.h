#ifndef WFTHELPER_H
#define WFTHELPER_H
#include "OneSiteSpaces.hh"
#include "Vector.h"

namespace Dmrg {

template <typename ModelType, typename VectorWithOffsetType, typename WaveFunctionTransfType>
class WftHelper {

public:

	typedef typename PsimagLite::Vector<VectorWithOffsetType*>::Type VectorVectorWithOffsetType;
	typedef typename ModelType::ModelHelperType                      ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType             LeftRightSuperType;
	typedef PsimagLite::Vector<SizeType>::Type                       VectorSizeType;
	using OneSiteSpacesType = OneSiteSpaces<ModelType>;

	WftHelper(const ModelType&              model,
	          const LeftRightSuperType&     lrs,
	          const WaveFunctionTransfType& wft)
	    : model_(model)
	    , lrs_(lrs)
	    , wft_(wft)
	{ }

	void
	wftSome(VectorVectorWithOffsetType& tvs, SizeType site, SizeType begin, SizeType end) const
	{
		for (SizeType index = begin; index < end; ++index) {
			assert(index < tvs.size());
			const VectorWithOffsetType& src = *tvs[index];
			if (src.size() == 0)
				continue;
			VectorWithOffsetType phiNew;
			wftOneVector(phiNew, src, site);
			*tvs[index] = phiNew;
		}
	}

	void wftOneVector(VectorWithOffsetType&       phiNew,
	                  const VectorWithOffsetType& src,
	                  SizeType                    site) const
	{
		phiNew.populateFromQns(src, lrs_.super());

		// OK, now that we got the partition number right, let's wft:
		ProgramGlobals::DirectionEnum dir
		    = ProgramGlobals::DirectionEnum::EXPAND_SYSTEM; // FIXME!
		OneSiteSpacesType oneSiteSpaces(site, dir, model_);
		wft_.setInitialVector(phiNew, src, lrs_, oneSiteSpaces);
	}

private:

	const ModelType&              model_;
	const LeftRightSuperType&     lrs_;
	const WaveFunctionTransfType& wft_;
};
}
#endif // WFTHELPER_H
