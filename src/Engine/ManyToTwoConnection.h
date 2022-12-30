#ifndef MANYTOTWOCONNECTION_H
#define MANYTOTWOCONNECTION_H
#include "ProgramGlobals.h"
#include "MetaOpForConnection.hh"

namespace Dmrg {

template<typename ModelLinksType,
         typename LeftRightSuperType,
         typename SuperOpHelperType>
class ManyToTwoConnection {

public:

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef std::pair<char, char> PairCharType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelLinksType::TermType::OneLinkType ModelTermLinkType;
	typedef typename ModelLinksType::HermitianEnum HermitianEnum;
	using PairMetaOpForConnection = std::pair<MetaOpForConnection, MetaOpForConnection>;

	ManyToTwoConnection(const VectorSizeType& hItems,
	                    ProgramGlobals::ConnectionEnum type,
	                    const ModelTermLinkType& oneLink,
	                    const LeftRightSuperType& lrs,
	                    const SuperOpHelperType& superOpHelper)
	    : oneLink_(oneLink), lrs_(lrs)
	{
		if (hItems.size() == 2) {
			assert(oneLink.indices.size() == 2);
			finalIndices_.first = locationFirst(hItems[0], oneLink.indices[0], type);
			finalIndices_.second = locationSecond(hItems[1], oneLink.indices[1], type);
			assert(oneLink.mods.size() == 2);
			mods_ = PairCharType(oneLink.mods[0], oneLink.mods[1]);
		} else {
			PairMetaOpForConnection finals = superOpHelper.finalIndices(hItems, type);
			mods_ = PairCharType('N', 'N');
			convertNonLocals(finals, type);
		}
	}

	const PairSizeType& finalIndices() const { return finalIndices_; }

	const PairCharType& finalMods() const { return mods_; }

	// WARNING: It doesn't consider the value of connection, so basically
	// it's OK only if value of connection is real
	bool connectionIsHermitian(const ModelLinksType& modelLinks) const
	{
		return (oneLink_.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION) ?
		            linkIsHermitianFermion(modelLinks) : linkIsHermitianBoson(modelLinks);
	}

private:

	bool linkIsHermitianFermion(const ModelLinksType& modelLinks) const
	{
		assert(oneLink_.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION);

		assert(oneLink_.indices.size() == 2);
		HermitianEnum h1 = modelLinks.getHermitianProperty(oneLink_.indices[0]);
		HermitianEnum h2 = modelLinks.getHermitianProperty(oneLink_.indices[1]);

		bool isHermit1 = (h1 == ModelLinksType::HERMIT_PLUS);
		bool isHermit2 = (h2 == ModelLinksType::HERMIT_PLUS);
		bool isAnti1 = (h1 == ModelLinksType::HERMIT_MINUS);
		bool isAnti2 = (h2 == ModelLinksType::HERMIT_MINUS);
		bool b1 = (isHermit1 && isAnti2);
		bool b2 = (isAnti1 && isHermit2);
		return (b1 || b2);
	}

	bool linkIsHermitianBoson(const ModelLinksType& modelLinks) const
	{
		assert(oneLink_.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::BOSON);

		assert(oneLink_.indices.size() == 2);
		HermitianEnum h1 = modelLinks.getHermitianProperty(oneLink_.indices[0]);
		HermitianEnum h2 = modelLinks.getHermitianProperty(oneLink_.indices[1]);

		bool isHermit1 = (h1 == ModelLinksType::HERMIT_PLUS);
		bool isHermit2 = (h2 == ModelLinksType::HERMIT_PLUS);

		return (isHermit1 && isHermit2);
	}

	SizeType locationFirst(SizeType hItems0,
	                       SizeType sigma,
	                       ProgramGlobals::ConnectionEnum type) const
	{
		const ProgramGlobals::SysOrEnvEnum sysOrEnv =
		        (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::SYSTEM : ProgramGlobals::SysOrEnvEnum::ENVIRON;

		SizeType i = PsimagLite::indexOrMinusOne(lrs_.super().block(), hItems0);

		const SizeType offset = lrs_.left().block().size();

		SizeType site1Corrected = (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            i : i - offset;

		return finalIndex(sysOrEnv, site1Corrected, sigma);
	}

	SizeType locationSecond(SizeType hItems1,
	                       SizeType sigma,
	                       ProgramGlobals::ConnectionEnum type) const
	{
		const ProgramGlobals::SysOrEnvEnum envOrSys =
		        (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::ENVIRON : ProgramGlobals::SysOrEnvEnum::SYSTEM;

		SizeType j = PsimagLite::indexOrMinusOne(lrs_.super().block(), hItems1);

		const SizeType offset = lrs_.left().block().size();

		SizeType site2Corrected = (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            j - offset : j;

		return finalIndex(envOrSys, site2Corrected, sigma);
	}

	SizeType finalIndex(ProgramGlobals::SysOrEnvEnum type,
	                    SizeType i,
	                    SizeType sigma) const
	{
		assert(type == ProgramGlobals::SysOrEnvEnum::SYSTEM ||
		       type == ProgramGlobals::SysOrEnvEnum::ENVIRON);

		return (type == ProgramGlobals::SysOrEnvEnum::SYSTEM) ?
		            lrs_.left(). localOperatorIndex(i, sigma) :
		            lrs_.right().localOperatorIndex(i, sigma);
	}

	void convertNonLocals(const PairMetaOpForConnection& pairMetas,
	                      ProgramGlobals::ConnectionEnum type)
	{
		if (pairMetas.first.site >= 0) {
			finalIndices_.first = locationFirst(pairMetas.first.site, pairMetas.first.index, type);
		} else {
			finalIndices_.first = pairMetas.first.index;
		}

		if (pairMetas.second.site >= 0) {
			finalIndices_.second = locationSecond(pairMetas.second.site, pairMetas.second.index, type);
		} else {
			finalIndices_.second = pairMetas.second.index;
		}
	}

	const ModelTermLinkType& oneLink_;
	const LeftRightSuperType& lrs_;
	PairSizeType finalIndices_;
	PairCharType mods_;
};
}
#endif // MANYTOTWOCONNECTION_H
