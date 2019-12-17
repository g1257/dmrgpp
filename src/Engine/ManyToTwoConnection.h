#ifndef MANYTOTWOCONNECTION_H
#define MANYTOTWOCONNECTION_H
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelTermLinkType, typename LeftRightSuperType>
class ManyToTwoConnection {

public:

	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef std::pair<char, char> PairCharType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	ManyToTwoConnection(const VectorSizeType& hItems,
	                    ProgramGlobals::ConnectionEnum type,
	                    const ModelTermLinkType& oneLink,
	                    const LeftRightSuperType& lrs)
	    : oneLink_(oneLink), lrs_(lrs)
	{
		finalIndices_ = finalIndices(hItems, type);
		assert(oneLink.mods.size() == 2);
		mods_ = PairCharType(oneLink.mods[0], oneLink.mods[1]);
	}

	const PairSizeType& finalIndices() const { return finalIndices_; }

	const PairCharType& finalMods() const { return mods_; }

private:


	PairSizeType finalIndices(const VectorSizeType& hItems,
	                          ProgramGlobals::ConnectionEnum type) const
	{
		const ProgramGlobals::SysOrEnvEnum sysOrEnv =
		        (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::SYSTEM : ProgramGlobals::SysOrEnvEnum::ENVIRON;
		const ProgramGlobals::SysOrEnvEnum envOrSys =
		        (type == ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON) ?
		            ProgramGlobals::SysOrEnvEnum::ENVIRON : ProgramGlobals::SysOrEnvEnum::SYSTEM;

		assert(hItems.size() == 2);
		assert(oneLink_.indices.size() == 2);
		PairSizeType finalIndex0;
		finalIndex0.first = finalIndex(sysOrEnv, hItems[0], oneLink_.indices[0]);
		finalIndex0.second = finalIndex(envOrSys, hItems[1], oneLink_.indices[1]);
		return finalIndex0;
	}

	SizeType finalIndex(ProgramGlobals::SysOrEnvEnum type,
	                    SizeType i,
	                    SizeType sigma) const
	{
		PairSizeType ii;
		if (type == ProgramGlobals::SysOrEnvEnum::SYSTEM) {
			ii = lrs_.left().getOperatorIndices(i, sigma);
		} else {
			assert(type == ProgramGlobals::SysOrEnvEnum::ENVIRON);
			ii = lrs_.right().getOperatorIndices(i, sigma);
		}

		return ii.first;
	}

	const ModelTermLinkType& oneLink_;
	const LeftRightSuperType& lrs_;
	PairSizeType finalIndices_;
	PairCharType mods_;
};
}
#endif // MANYTOTWOCONNECTION_H
