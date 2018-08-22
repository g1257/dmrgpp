/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
/** \ingroup DMRG */
/*@{*/

/*! \file DmrgWaveStruct.h
 *
 *  DOC NEEDED FIXME (This file should go in Wft/ directory perhaps)
 */
#ifndef DMRG_WAVE_H
#define DMRG_WAVE_H
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg {

template<typename LeftRightSuperType_>
struct DmrgWaveStruct {

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	DmrgWaveStruct()
	    : lrs_("pSE", "pSprime", "pEprime") {}

	void setLrs(const LeftRightSuperType& lrs)
	{
		lrs_.dontCopyOperators(lrs);
	}

	void setTransform(const BlockDiagonalMatrixType& m, SizeType what)
	{
		if (what == ProgramGlobals::SYSTEM) {
			ws_ = m;
			return;
		}

		assert(what == ProgramGlobals::ENVIRON);
		we_ = m;
	}

	const BlockDiagonalMatrixType& getTransform(SizeType what) const
	{
		if (what == ProgramGlobals::SYSTEM)
			return ws_;

		assert(what == ProgramGlobals::ENVIRON);
		return we_;
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

	template<typename IoInputType>
	void read(IoInputType& io,
	          PsimagLite::String prefix,
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsInputLike<IoInputType>::True, int>::Type = 0)
	{
		io.read(ws_, prefix + "/Ws");
		io.read(we_, prefix + "/We");
		lrs_.read(io, prefix);
	}

	template<typename IoOutputType>
	void write(IoOutputType& io,
	           PsimagLite::String prefix,
	           typename PsimagLite::EnableIf<
	           PsimagLite::IsOutputLike<IoOutputType>::True, int>::Type = 0) const
	{
		io.createGroup(prefix);
		io.write(ws_, prefix + "/Ws");
		io.write(we_, prefix + "/We");
		lrs_.write(io, prefix, LeftRightSuperType::SAVE_ALL, false);
	}

private:

	BlockDiagonalMatrixType ws_;
	BlockDiagonalMatrixType we_;
	LeftRightSuperType lrs_;
}; // struct DmrgWaveStruct

} // namespace Dmrg 

/*@}*/
#endif // DMRG_WAVE_H
