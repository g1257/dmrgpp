/*
Copyright (c) 2009,-2012 UT-Battelle, LLC
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
/** \file ParallelDensityMatrix.h
*/

#ifndef PARALLEL_DENSITY_MATRIX_H
#define PARALLEL_DENSITY_MATRIX_H

#include "ProgramGlobals.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename BlockMatrixType,typename BasisWithOperatorsType,
         typename TargettingType>
class ParallelDensityMatrix {

	typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;
	typedef typename TargettingType::VectorWithOffsetType TargetVectorType;
	typedef typename TargetVectorType::value_type DensityMatrixElementType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef PsimagLite::Concurrency ConcurrencyType;

public:

	typedef typename PsimagLite::Real<DensityMatrixElementType>::Type RealType;

	ParallelDensityMatrix(const TargettingType& target,
	                      const BasisWithOperatorsType& pBasis,
	                      const BasisWithOperatorsType& pBasisSummed,
	                      const BasisType& pSE,
	                      int direction,
	                      SizeType m,
	                      BuildingBlockType& matrixBlock)
	    : target_(target),
	      pBasis_(pBasis),
	      pBasisSummed_(pBasisSummed),
	      pSE_(pSE),
	      direction_(direction),
	      m_(m),
	      matrixBlock_(matrixBlock),
	      hasMpi_(PsimagLite::Concurrency::hasMpi())
	{}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      pthread_mutex_t*)
	{
		SizeType mpiRank = (hasMpi_) ?
		            PsimagLite::MPI::commRank(PsimagLite::MPI::COMM_WORLD) : 0;
		SizeType npthreads = PsimagLite::Concurrency::npthreads;

		ConcurrencyType::mpiDisableIfNeeded(mpiRank,
		                                    blockSize,
		                                    "ParallelDensityMatrix",
		                                    total);

		for (SizeType p=0;p<blockSize;p++) {
			SizeType ix = (threadNum+npthreads*mpiRank)*blockSize + p;
			if (ix>=target_.size()) break;
			RealType wnorm = target_.normSquared(ix);
			if (fabs(wnorm) < 1e-6) continue;
			RealType w = target_.weight(ix)/wnorm;
			initPartition(matrixBlock_,pBasis_,m_,target_(ix),
			              pBasisSummed_,pSE_,direction_,w);
		}
	}

	void initPartition(BuildingBlockType& matrixBlock,
	                   BasisWithOperatorsType const &pBasis,
	                   SizeType m,
	                   const TargetVectorType& v,
	                   BasisWithOperatorsType const &pBasisSummed,
	                   BasisType const &pSE,
	                   SizeType direction,
	                   RealType weight)
	{
		if (direction!=ProgramGlobals::EXPAND_SYSTEM)
			initPartitionExpandEnviron(matrixBlock,pBasis,m,v,pBasisSummed,pSE,weight);
		else
			initPartitionExpandSystem(matrixBlock,pBasis,m,v,pBasisSummed,pSE,weight);
	}

private:

	void initPartitionExpandEnviron(BuildingBlockType& matrixBlock,
	                                BasisWithOperatorsType const &pBasis,
	                                SizeType m,
	                                const TargetVectorType& v,
	                                BasisWithOperatorsType const &pBasisSummed,
	                                BasisType const &pSE,
	                                RealType weight)
	{

		SizeType ns=pBasisSummed.size();
		SizeType ne=pSE.size()/ns;
		SizeType start = pBasis.partition(m);
		SizeType length = pBasis.partition(m+1) - start;
		SizeType total=pBasisSummed.size();

		for (SizeType i=0;i<length;++i) {
			SizeType ieff = i +start;
			for (SizeType j=0;j<length;++j) {
				matrixBlock(i,j) += densityMatrixExpandEnviron(ieff,
				                                               j+start,
				                                               v,
				                                               total,
				                                               pSE,
				                                               ns,
				                                               ne)*weight;
			}
		}
	}

	void initPartitionExpandSystem(BuildingBlockType& matrixBlock,
	                               BasisWithOperatorsType const &pBasis,
	                               SizeType m,
	                               const TargetVectorType& v,
	                               BasisWithOperatorsType const &pBasisSummed,
	                               BasisType const &pSE,
	                               RealType weight)
	{
		SizeType ne = pBasisSummed.size();
		SizeType ns = pSE.size()/ne;
		SizeType start = pBasis.partition(m);
		SizeType length = pBasis.partition(m+1) - start;
		SizeType total=pBasisSummed.size();

		for (SizeType i=0;i<length;++i) {
			SizeType ieff = i +start;
			for (SizeType j=0;j<length;++j) {
				matrixBlock(i,j) += densityMatrixExpandSystem(ieff,
				                                              j+start,
				                                              v,
				                                              total,
				                                              pSE,
				                                              ns,
				                                              ne)*weight;
			}
		}
	}

	DensityMatrixElementType densityMatrixExpandEnviron(
	        SizeType alpha1,
	        SizeType alpha2,
	        const TargetVectorType& v,
	        SizeType total,
	        BasisType const &pSE,
	        SizeType ns,
	        SizeType)
	{
		DensityMatrixElementType sum=0;

		SizeType x2 = alpha2*ns;
		SizeType x1 = alpha1*ns;
		for (SizeType beta=0;beta<total;beta++) {
			SizeType ii = pSE.permutationInverse(beta + x1);
			int sector1 = v.index2Sector(ii);
			if (sector1 < 0) continue;
			SizeType start1 = v.offset(sector1);

			SizeType jj = pSE.permutationInverse(beta + x2);
			int sector2 = v.index2Sector(jj);
			if (sector2 < 0) continue;
			SizeType start2 = v.offset(sector2);

			sum += v.fastAccess(sector1,ii-start1)*
			        std::conj(v.fastAccess(sector2,jj-start2));
		}
		return sum;
	}

	DensityMatrixElementType densityMatrixExpandSystem(
	        SizeType alpha1,
	        SizeType alpha2,
	        const TargetVectorType& v,
	        SizeType total,
	        BasisType const &pSE,
	        SizeType ns,
	        SizeType)
	{
		DensityMatrixElementType sum=0;

		SizeType totalNs = total * ns;

		for (SizeType betaNs=0;betaNs<totalNs;betaNs+=ns) {
			SizeType ii = pSE.permutationInverse(alpha1+betaNs);
			int sector1 = v.index2Sector(ii);
			if (sector1 < 0) continue;
			SizeType start1 = v.offset(sector1);

			SizeType jj = pSE.permutationInverse(alpha2+betaNs);
			int sector2 = v.index2Sector(jj);
			if (sector2 < 0) continue;
			SizeType start2 = v.offset(sector2);

			sum += v.fastAccess(sector1,ii-start1)*
			        std::conj(v.fastAccess(sector2,jj-start2));
		}

		return sum;
	}

	const TargettingType& target_;
	const BasisWithOperatorsType& pBasis_;
	const BasisWithOperatorsType& pBasisSummed_;
	const BasisType& pSE_;
	int direction_;
	SizeType m_;
	BuildingBlockType& matrixBlock_;
	bool hasMpi_;
}; // class ParallelDensityMatrix
} // namespace Dmrg

/*@}*/
#endif // PARALLEL_DENSITY_MATRIX_H

