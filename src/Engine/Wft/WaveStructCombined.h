#ifndef WAVESTRUCTCOMBINED_H
#define WAVESTRUCTCOMBINED_H
#include "Io/IoNg.h"
#include "WaveStructSvd.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename LeftRightSuperType_>
class WaveStructCombined {

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
	typedef typename PsimagLite::Stack<WaveStructSvdType>::Type WftStackType;

	WaveStructCombined()
	    : lrs_("pSE", "pSprime", "pEprime"), needsPop_(false)
	{}

	void read(PsimagLite::IoNg::In& io, PsimagLite::String prefix)
	{
		lrs_.read(io, prefix);
		io.read(wsStack_, prefix + "/wsStack");
		io.read(weStack_, prefix + "/weStack");
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix) const
	{
		writePartial(io, prefix);
		WftStackType wsStack = wsStack_;
		io.write(wsStack, prefix + "/wsStack");
		WftStackType weStack = weStack_;
		io.write(weStack, prefix + "/weStack");
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix)
	{
		writePartial(io, prefix);
		io.write(wsStack_, prefix + "/wsStack");
		io.write(weStack_, prefix + "/weStack");
	}

	void beforeWft(ProgramGlobals::DirectionEnum dir,
	               bool twoSiteDmrg,
	               SizeType counter)
	{
		WftStackType& stack = (dir == ProgramGlobals::EXPAND_ENVIRON) ? wsStack_ :
		                                                                weStack_;
		const PsimagLite::String label = (dir == ProgramGlobals::EXPAND_ENVIRON) ?
		            "system" : "environ";

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
		if (stack.size() > 0 && counter == 0)
			stack.pop();
		else
			needsPop_ = true;
	}

	void afterWft(ProgramGlobals::DirectionEnum dir)
	{
		WftStackType& stack = (dir == ProgramGlobals::EXPAND_ENVIRON) ? wsStack_ :
		                                                                weStack_;
		const PsimagLite::String label = (dir == ProgramGlobals::EXPAND_ENVIRON) ?
		            "system" : "environ";

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
		case ProgramGlobals::INFINITE:
			if (direction == ProgramGlobals::EXPAND_SYSTEM) {
				wsStack_.push(wave);
			} else {
				weStack_.push(wave);
			}

			break;
		case ProgramGlobals::EXPAND_ENVIRON:
			if (direction != ProgramGlobals::EXPAND_ENVIRON)
				err("EXPAND_ENVIRON but option==0\n");
			weStack_.push(wave);
			break;
		case ProgramGlobals::EXPAND_SYSTEM:
			if (direction != ProgramGlobals::EXPAND_SYSTEM)
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
		assert(sysOrEnv == ProgramGlobals::SYSTEM || weStack_.size() > 0);
		assert(sysOrEnv != ProgramGlobals::SYSTEM || wsStack_.size() > 0);
		return (sysOrEnv == ProgramGlobals::SYSTEM) ? wsStack_.top() : weStack_.top();
	}

	const LeftRightSuperType& lrs() const
	{
		return lrs_;
	}

	const BlockDiagonalMatrixType& getTransform(ProgramGlobals::SysOrEnvEnum dir) const
	{
		return getWave(dir).u();
	}

private:

	void writePartial(PsimagLite::IoSelector::Out& io, PsimagLite::String prefix) const
	{
		io.createGroup(prefix);
		lrs_.write(io, prefix, LeftRightSuperType::SAVE_ALL, false);
	}

	LeftRightSuperType lrs_;
	WftStackType wsStack_;
	WftStackType weStack_;
	bool needsPop_;
};
}
#endif // WAVESTRUCTCOMBINED_H
