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

	void read(PsimagLite::IoNg::In& io, PsimagLite::String prefix2)
	{
		PsimagLite::String prefix = prefix2 + "/WaveStructCombined";
		lrs_.read(io, prefix);
		left_.read(io, prefix);
		right_.read(io, prefix);
	}

	void write(PsimagLite::IoNg::Out& io, PsimagLite::String prefix2) const
	{
		io.createGroup(prefix2);
		PsimagLite::String prefix = prefix2 + "/WaveStructCombined";
		lrs_.write(io, prefix, LeftRightSuperType::SAVE_ALL, false);
		left_.write(io, prefix);
		right_.write(io, prefix);
	}

	void setLrs(const LeftRightSuperType& lrs)
	{
		lrs_.dontCopyOperators(lrs);
	}

	void clear()
	{
		left_.clear();
		right_.clear();
	}

	void setWave(const WaveStructSvdType& wave,
	             ProgramGlobals::SysOrEnvEnum sysOrEnv)
	{
		(sysOrEnv == ProgramGlobals::SYSTEM) ? left_ : right_ = wave;
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
