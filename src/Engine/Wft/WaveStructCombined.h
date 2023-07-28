#ifndef WAVESTRUCTCOMBINED_H
#define WAVESTRUCTCOMBINED_H
#include "../DiskOrMemoryStack.h"
#include "BasisTraits.hh"
#include "Io/IoNg.h"
#include "ProgramGlobals.h"
#include "WaveStructSvd.h"

namespace Dmrg
{

template <typename LeftRightSuperType_>
class WaveStructCombined
{

public:

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef WaveStructSvd<LeftRightSuperType> WaveStructSvdType;
	typedef typename WaveStructSvdType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename WaveStructSvdType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename WaveStructSvdType::VectorVectorRealType VectorVectorRealType;
	typedef typename WaveStructSvdType::VectorMatrixType VectorMatrixType;
	typedef typename WaveStructSvdType::VectorQnType VectorQnType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisType::BlockType VectorSizeType;
	typedef DiskOrMemoryStack<WaveStructSvdType> WftStackType;

	WaveStructCombined(bool onDisk,
	    const PsimagLite::String filename,
	    const BasisTraits& basisTraits)
	    : lrs_("pSE", "pSprime", "pEprime", basisTraits)
	    , wsStack_(onDisk, filename, "Wstacks", "system", basisTraits)
	    , weStack_(onDisk, filename, "Wstacks", "environ", basisTraits)
	    , needsPop_(false)
	{
	}

	void read(PsimagLite::IoNg::In& io, PsimagLite::String prefix)
	{
		lrs_.read(io, prefix);
		io.read(wsStack_, prefix + "/wsStack");
		io.read(weStack_, prefix + "/weStack");
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix) const
	{
		writePartial(io, prefix);
		wsStack_.write(prefix + "/wsStack", io.serializer());
		weStack_.write(prefix + "/weStack", io.serializer());
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix)
	{
		writePartial(io, prefix);
		wsStack_.write(prefix + "/wsStack", io.serializer());
		weStack_.write(prefix + "/weStack", io.serializer());
	}

	void beforeWft(ProgramGlobals::DirectionEnum dir,
	    bool twoSiteDmrg,
	    bool bounce)
	{
		WftStackType& stack = (dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) ? wsStack_
											     : weStack_;
		const PsimagLite::String label = (dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) ? "system" : "environ";

		needsPop_ = false;

		if (twoSiteDmrg) {
			if (stack.size() == 0)
				err("Stack for " + label + " is empty\n");
			if (stack.size() > 1)
				stack.pop();
			else
				needsPop_ = true;
			return;
		}

		assert(!twoSiteDmrg);
		if (stack.size() > 1 && bounce)
			stack.pop();
		else
			needsPop_ = true;
	}

	void afterWft(ProgramGlobals::DirectionEnum dir)
	{
		WftStackType& stack = (dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) ? wsStack_
											     : weStack_;
		const PsimagLite::String label = (dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) ? "system" : "environ";

		if (!needsPop_)
			return;

		assert(needsPop_);
		if (stack.size() == 0)
			err("Stack for " + label + " is empty\n");

		stack.pop();
		needsPop_ = false;
	}

	void push(const BlockDiagonalMatrixType& transform,
	    ProgramGlobals::DirectionEnum direction,
	    const VectorMatrixType& vts,
	    const VectorVectorRealType& s,
	    const VectorQnType& qns,
	    ProgramGlobals::DirectionEnum dir)
	{
		WaveStructSvdType wave(transform, vts, s, qns);

		switch (dir) {
		case ProgramGlobals::DirectionEnum::INFINITE:
			if (direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				wsStack_.push(wave);
			} else {
				weStack_.push(wave);
			}

			break;
		case ProgramGlobals::DirectionEnum::EXPAND_ENVIRON:
			if (direction != ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
				err("EXPAND_ENVIRON but option==0\n");
			weStack_.push(wave);
			break;
		case ProgramGlobals::DirectionEnum::EXPAND_SYSTEM:
			if (direction != ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
				err("EXPAND_SYSTEM but option==1\n");
			wsStack_.push(wave);
			break;
		}
	}

	void setLrs(const LeftRightSuperType& lrs)
	{
		lrs_.dontCopyOperators(lrs);
	}

	const WaveStructSvdType& getWave(ProgramGlobals::SysOrEnvEnum sysOrEnv) const
	{
		assert(sysOrEnv == ProgramGlobals::SysOrEnvEnum::SYSTEM || weStack_.size() > 0);
		assert(sysOrEnv != ProgramGlobals::SysOrEnvEnum::SYSTEM || wsStack_.size() > 0);
		return (sysOrEnv == ProgramGlobals::SysOrEnvEnum::SYSTEM) ? wsStack_.top()
									  : weStack_.top();
	}

	const LeftRightSuperType& lrs() const
	{
		return lrs_;
	}

	const BlockDiagonalMatrixType& getTransform(ProgramGlobals::SysOrEnvEnum dir) const
	{
		return getWave(dir).u();
	}

	SizeType size(ProgramGlobals::SysOrEnvEnum sysOrEnv) const
	{
		return (sysOrEnv == ProgramGlobals::SysOrEnvEnum::SYSTEM) ? wsStack_.size()
									  : weStack_.size();
	}

	const BlockDiagonalMatrixType& multiPointGetTransform(SizeType ind,
	    ProgramGlobals::DirectionEnum dir) const
	{
		return (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? wsStack_[ind].u()
									     : weStack_[ind].u();
	}

private:

	void writePartial(PsimagLite::IoSelector::Out& io, PsimagLite::String prefix) const
	{
		io.createGroup(prefix);
		lrs_.write(io, prefix, BasisWithOperatorsType::SaveEnum::ALL, false);
	}

	LeftRightSuperType lrs_;
	WftStackType wsStack_;
	WftStackType weStack_;
	bool needsPop_;
};
}
#endif // WAVESTRUCTCOMBINED_H
