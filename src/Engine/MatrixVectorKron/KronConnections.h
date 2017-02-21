/*
Copyright (c) 2012, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file KronConnections.h
 *
 *
 */

#ifndef KRON_CONNECTIONS_H
#define KRON_CONNECTIONS_H

#include "Matrix.h"
#include "Concurrency.h"
#include "InitKron.h"

namespace Dmrg {

template<typename ModelType, typename ModelHelperType_>
class KronConnections {
	typedef InitKron<ModelType,ModelHelperType_> InitKronType;

	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename InitKronType::GenGroupType GenGroupType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename InitKronType::RealType RealType;


	KronConnections(const InitKronType& initKron,
	                VectorType& x,
	                const VectorType& y,
	                SizeType nthreads)
	    : initKron_(initKron),x_(x),y_(y)
	{}

	SizeType tasks() const
	{
		assert(initKron_.patch(GenIjPatchType::LEFT).size() ==
		       initKron_.patch(GenIjPatchType::RIGHT).size());
		return initKron_.patch(GenIjPatchType::LEFT).size();
	}

	void doTask(SizeType outPatch, SizeType)
	{
		SizeType nC = initKron_.connections();
		SizeType total = tasks();
		for (SizeType inPatch=0;inPatch<total;++inPatch) {
			for (SizeType ic=0;ic<nC;++ic) {
				const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
				const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

				SizeType i = initKron_.patch(GenIjPatchType::LEFT)[inPatch];
				SizeType j = initKron_.patch(GenIjPatchType::RIGHT)[inPatch];

				SizeType ip = initKron_.patch(GenIjPatchType::LEFT)[outPatch];
				SizeType jp = initKron_.patch(GenIjPatchType::RIGHT)[outPatch];


				const MatrixDenseOrSparseType& tmp1 =  xiStruct(ip,i);
				const MatrixDenseOrSparseType& tmp2 =  yiStruct(j,jp);

				kronMult(x_, &(y_[0]), 'n', 'n', tmp1, tmp2);
			}
		}

	}

	void sync()
	{}

	const InitKronType& initKron_;
	VectorType& x_;
	const VectorType& y_;
}; //class KronConnections

} // namespace PsimagLite

/*@}*/

#endif // KRON_CONNECTIONS_H
