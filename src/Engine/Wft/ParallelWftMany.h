/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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
/** \file ParallelWftMany.h
 */

#ifndef DMRG_PARALLEL_WFT_MANY_H
#define DMRG_PARALLEL_WFT_MANY_H

#include "Concurrency.h"
#include "Vector.h"

namespace Dmrg
{

template <typename VectorWithOffsetType,
    typename WaveFunctionTransfType,
    typename LeftRightSuperType>
class ParallelWftMany
{

	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type
	    VectorVectorWithOffsetType;

public:

	typedef typename VectorWithOffsetType::value_type VectorElementType;
	typedef typename PsimagLite::Real<VectorElementType>::Type RealType;

	ParallelWftMany(VectorVectorWithOffsetType& targetVectors,
	    SizeType nk,
	    const WaveFunctionTransfType& wft,
	    const LeftRightSuperType& lrs)
	    : targetVectors_(targetVectors)
	    , nk_(nk)
	    , wft_(wft)
	    , lrs_(lrs)
	{
	}

	void thread_function_(SizeType threadNum,
	    SizeType blockSize,
	    SizeType total,
	    ConcurrencyType::MutexType*)
	{
		SizeType nk = nk_;
		SizeType mpiRank = PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD);
		SizeType npthreads = PsimagLite::Concurrency::npthreads;

		ConcurrencyType::mpiDisableIfNeeded(mpiRank, blockSize, "ParallelWftMany", total);

		for (SizeType p = 0; p < blockSize; p++) {
			SizeType ix = (threadNum + npthreads * mpiRank) * blockSize + p + 1;
			if (ix >= targetVectors_.size())
				break;
			VectorWithOffsetType phiNew = targetVectors_[0];
			wft_.setInitialVector(phiNew, targetVectors_[ix], lrs_, nk);
			targetVectors_[ix] = phiNew;
		}
	}

private:

	VectorVectorWithOffsetType& targetVectors_;
	SizeType nk_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
}; // class ParallelWftMany
} // namespace Dmrg

/*@}*/
#endif // DMRG_PARALLEL_WFT_MANY_H
