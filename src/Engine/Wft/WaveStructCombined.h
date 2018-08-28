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

	WaveStructCombined()
	    : lrs_("pSE", "pSprime", "pEprime")
	{}

	void read(PsimagLite::IoNg::In& io, PsimagLite::String prefix)
	{
		lrs_.read(io, prefix);
		left_.read(io, prefix + "/left");
		right_.read(io, prefix + "/right");
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix) const
	{
		io.createGroup(prefix);
		lrs_.write(io, prefix, LeftRightSuperType::SAVE_ALL, false);
		left_.write(io, prefix + "/left");
		right_.write(io, prefix + "/right");
	}

	void setLrs(const LeftRightSuperType& lrs)
	{
		lrs_.dontCopyOperators(lrs);
	}

	void setWave(const WaveStructSvdType& wave,
	             ProgramGlobals::SysOrEnvEnum sysOrEnv)
	{
		if (sysOrEnv == ProgramGlobals::SYSTEM) left_  = wave;
		else right_ = wave;
	}

	const WaveStructSvdType& getWave(ProgramGlobals::SysOrEnvEnum sysOrEnv) const
	{
		return (sysOrEnv == ProgramGlobals::SYSTEM) ? left_ : right_;
	}

	const LeftRightSuperType& lrs() const
	{
		return lrs_;
	}

	const BlockDiagonalMatrixType& getTransform(SizeType dir) const
	{

		return (dir == ProgramGlobals::SYSTEM) ? left_.u() :
		                                         right_.u();
	}

private:

	LeftRightSuperType lrs_;
	WaveStructSvdType left_;
	WaveStructSvdType right_;
};
}
#endif // WAVESTRUCTCOMBINED_H
