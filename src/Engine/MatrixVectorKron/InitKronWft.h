/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file InitKronWft.h
 *
 *
 */
#ifndef INITKRON_WFT_H
#define INITKRON_WFT_H
#include "ProgramGlobals.h"
#include "InitKronBase.h"
#include "Vector.h"

namespace Dmrg {

template<typename LeftRightSuperType_, typename WftOptionsType, typename DmrgWaveStructType>
class InitKronWft : public InitKronBase<LeftRightSuperType_> {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef InitKronBase<LeftRightSuperType> BaseType;
	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::LinkType LinkType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BaseType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;
	typedef typename LinkType::PairSizeType PairSizetype;
	typedef typename LinkType::PairCharType PairCharType;

	InitKronWft(const LeftRightSuperType& lrs,
	            SizeType m,
	            SizeType qn,
	            const WftOptionsType& wftOptions,
	            const DmrgWaveStructType& dmrgWaveStruct)
	    : BaseType(lrs, m, qn, wftOptions.denseSparseThreshold),
	      wftOptions_(wftOptions)
	{
		// FIXME: dmrgWaveStruct_.lrs() vs. modelHelper_.lrs(), which one to use?
		SparseMatrixType we;
		dmrgWaveStruct.we.toSparse(we);
		const PairCharType nn('N', 'N');
		const PairSizetype dummy(0, 0);
		const ComplexOrRealType value = 1.0;

		LinkType link(0,
		              0,
		              ProgramGlobals::SYSTEM_ENVIRON,
		              value,
		              1,
		              ProgramGlobals::BOSON,
		              dummy,
		              nn,
		              1,
		              1,
		              0);

		SparseMatrixType ws;
		dmrgWaveStruct.ws.toSparse(ws);

		if (wftOptions_.dir == ProgramGlobals::EXPAND_SYSTEM) {
			this->addOneConnection(ws, we, link);
		} else {
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			SparseMatrixType wsT;
			transposeConjugate(wsT,ws);
			this->addOneConnection(wsT, weT, link);
		}
	}

	bool loadBalance() const
	{
		return wftOptions_.kronLoadBalance;
	}

private:

	const WftOptionsType& wftOptions_;

	InitKronWft(const InitKronWft&);

	InitKronWft& operator=(const InitKronWft&);
};
} // namespace Dmrg

#endif // INITKRON_WFT_H
