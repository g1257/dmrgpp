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

/*! \file InitKronHamiltonian.h
 *
 *
 */
#ifndef PREINITKRON_WFT_H
#define PREINITKRON_WFT_H
#include "ProgramGlobals.h"
#include "PreInitKronBase.h"
#include "Vector.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class PreInitKronWft : public PreInitKronBase<LeftRightSuperType> {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef PreInitKronBase<LeftRightSuperType> BaseType;
	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename LeftRightSuperType::LinkType LinkType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BaseType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;
	typedef typename LinkType::PairSizeType PairSizetype;
	typedef typename LinkType::PairCharType PairCharType;

	PreInitKronWft(const LeftRightSuperType& lrs)
	    : BaseType(lrs)
	{
		// FIXME: dmrgWaveStruct_.lrs() vs. modelHelper_.lrs(), which one to use?
		SparseMatrixType we;
		// dmrgWaveStruct_.we.toSparse(we);
		const PairCharType nn('N', 'N');
		const PairSizetype dummy(0, 0);
		const ComplexOrRealType value = 1.0;

		LinkType(0,
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
		//dmrgWaveStruct_.ws.toSparse(ws);
		LinkType link();
		//if (dir_ == ProgramGlobals::EXPAND_SYSTEM) {
		SparseMatrixType wsT;
		transposeConjugate(wsT,ws);
		addOneConnection(wsT, we, link);
		/*} else {
			SparseMatrixType weT;
			transposeConjugate(weT,we);
			addOneConnection(ws, weT, link);
		}*/
	}

private:

	PreInitKronWft(const PreInitKronWft&);

	PreInitKronWft& operator=(const PreInitKronWft&);
};
} // namespace Dmrg

#endif // PREINITKRON_WFT_H
