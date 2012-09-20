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

/*! \file ParallelDensityMatrix.h
 *
 * DOC TBW FIXME
 */
#ifndef PARALLEL_DENSITY_MATRIX_H
#define PARALLEL_DENSITY_MATRIX_H

#include "ProgramGlobals.h"


namespace Dmrg {
	
	template<typename RealType_,typename BlockMatrixType,typename BasisWithOperatorsType,
			 typename TargettingType>
	class ParallelDensityMatrix {

		typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;
		typedef typename TargettingType::VectorWithOffsetType TargetVectorType;
		typedef typename TargetVectorType::value_type DensityMatrixElementType;
		typedef typename BasisWithOperatorsType::BasisType BasisType;

	public:

		typedef RealType_ RealType;

		ParallelDensityMatrix(const TargettingType& target,
							  const BasisWithOperatorsType& pBasis,
							  const BasisWithOperatorsType& pBasisSummed,
							  const BasisType& pSE,
							  int direction,
							  size_t m,
							  BuildingBlockType& matrixBlock)
		: target_(target),
		  pBasis_(pBasis),
		  pBasisSummed_(pBasisSummed),
		  pSE_(pSE),
		  direction_(direction),
		  m_(m),
		  matrixBlock_(matrixBlock)
		{}

		void thread_function_(size_t threadNum,size_t blockSize,size_t total,pthread_mutex_t* myMutex)
		{
			for (size_t p=0;p<blockSize;p++) {
				size_t ix = threadNum * blockSize + p;
				if (ix>=target_.size()) break;
				RealType w = target_.weight(ix)/target_.normSquared(ix);
				initPartition(matrixBlock_,pBasis_,m_,target_(ix),
							  pBasisSummed_,pSE_,direction_,w);
			}
		}

		void initPartition(BuildingBlockType& matrixBlock,
						   BasisWithOperatorsType const &pBasis,
						   size_t m,
						   const TargetVectorType& v,
						   BasisWithOperatorsType const &pBasisSummed,
						   BasisType const &pSE,
						   size_t direction,
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
										size_t m,
										const TargetVectorType& v,
										BasisWithOperatorsType const &pBasisSummed,
										BasisType const &pSE,
										RealType weight)
		{

			size_t ns=pBasisSummed.size();
			size_t ne=pSE.size()/ns;

			for (size_t i=pBasis.partition(m);i<pBasis.partition(m+1);i++) {
				for (size_t j=pBasis.partition(m);j<pBasis.partition(m+1);j++) {

					matrixBlock(i-pBasis.partition(m),j-pBasis.partition(m)) +=
							densityMatrixExpandEnviron(i,j,v,pBasisSummed,pSE,ns,ne)*weight;
				}
			}
		}

		void initPartitionExpandSystem(BuildingBlockType& matrixBlock,
									   BasisWithOperatorsType const &pBasis,
									   size_t m,
									   const TargetVectorType& v,
									   BasisWithOperatorsType const &pBasisSummed,
									   BasisType const &pSE,
									   RealType weight)
		{
			size_t ne = pBasisSummed.size();
			size_t ns = pSE.size()/ne;

			for (size_t i=pBasis.partition(m);i<pBasis.partition(m+1);i++) {
				for (size_t j=pBasis.partition(m);j<pBasis.partition(m+1);j++) {

					matrixBlock(i-pBasis.partition(m),j-pBasis.partition(m)) +=
							densityMatrixExpandSystem(i,j,v,pBasisSummed,pSE,ns,ne)*weight;
				}
			}
		}

		DensityMatrixElementType densityMatrixExpandEnviron(
			size_t alpha1,
			size_t alpha2,
			const TargetVectorType& v,
			BasisWithOperatorsType const &pBasisSummed,
			BasisType const &pSE,
			size_t ns,
			size_t ne)
		{

			size_t total=pBasisSummed.size();

			DensityMatrixElementType sum=0;

			size_t x2 = alpha2*ns;
			size_t x1 = alpha1*ns;
			for (size_t beta=0;beta<total;beta++) {
				size_t jj = pSE.permutationInverse(beta + x2);
				size_t ii = pSE.permutationInverse(beta + x1);
				sum += v[ii] * std::conj(v[jj]);
			}
			return sum;
		}

		DensityMatrixElementType densityMatrixExpandSystem(
			size_t alpha1,
			size_t alpha2,
			const TargetVectorType& v,
			BasisWithOperatorsType const &pBasisSummed,
			BasisType const &pSE,
			size_t ns,
			size_t ne)
		{

			size_t total=pBasisSummed.size();

			DensityMatrixElementType sum=0;

			for (size_t beta=0;beta<total;beta++) {
				size_t jj = pSE.permutationInverse(alpha2+beta*ns);
				size_t ii = pSE.permutationInverse(alpha1+beta*ns);
				sum += v[ii] * std::conj(v[jj]);
			}
			return sum;
		}

		const TargettingType& target_;
		const BasisWithOperatorsType& pBasis_;
		const BasisWithOperatorsType& pBasisSummed_;
		const BasisType& pSE_;
		int direction_;
		size_t m_;
		BuildingBlockType& matrixBlock_;
	}; // class ParallelDensityMatrix
} // namespace Dmrg 

/*@}*/
#endif // PARALLEL_DENSITY_MATRIX_H
