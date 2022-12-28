#ifndef MANYTOTWOCONNECTION_H
#define MANYTOTWOCONNECTION_H
#include "ProgramGlobals.h"

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

	ManyToTwoConnection(const VectorSizeType& hItems,
	                    ProgramGlobals::ConnectionEnum type,
	                    const ModelTermLinkType& oneLink,
	                    const LeftRightSuperType& lrs,
	                    const SuperOpHelperType& superOpHelper)
	    : oneLink_(oneLink), lrs_(lrs)
	{
		finalIndices_ = finalIndices(hItems, type, superOpHelper);
		assert(oneLink.mods.size() == 2);
		mods_ = PairCharType(oneLink.mods[0], oneLink.mods[1]);
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

	PairSizeType finalIndices(const VectorSizeType& hItems,
	                          ProgramGlobals::ConnectionEnum type,
	                          const SuperOpHelperType& superOpHelper) const
	{
		if (hItems.size() == 4) {
			assert(type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON);
			return superOpHelper.finalIndices4sites(hItems, type);
		}

		const ProgramGlobals::SysOrEnvEnum sysOrEnv =
		        (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::SYSTEM : ProgramGlobals::SysOrEnvEnum::ENVIRON;
		const ProgramGlobals::SysOrEnvEnum envOrSys =
		        (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::ENVIRON : ProgramGlobals::SysOrEnvEnum::SYSTEM;

		assert(hItems.size() == 2);

		SizeType i = PsimagLite::indexOrMinusOne(lrs_.super().block(), hItems[0]);
		SizeType j = PsimagLite::indexOrMinusOne(lrs_.super().block(), hItems[1]);

		const SizeType offset = lrs_.left().block().size();

		SizeType site1Corrected = (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            i : i - offset;
		SizeType site2Corrected = (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            j - offset : j;

		assert(oneLink_.indices.size() == 2);
		PairSizeType finalIndex0;
		finalIndex0.first = finalIndex(sysOrEnv, site1Corrected, oneLink_.indices[0]);
		finalIndex0.second = finalIndex(envOrSys, site2Corrected, oneLink_.indices[1]);
		return finalIndex0;
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

	const ModelTermLinkType& oneLink_;
	const LeftRightSuperType& lrs_;
	PairSizeType finalIndices_;
	PairCharType mods_;
};
}
#endif // MANYTOTWOCONNECTION_H
