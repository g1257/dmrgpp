
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

#ifndef DENSITY_MATRIX_LOCAL_H
#define DENSITY_MATRIX_LOCAL_H
#include "ProgressIndicator.h"
#include "TypeToString.h"
#include "BlockDiagonalMatrix.h"
#include "DensityMatrixBase.h"
#include "ParallelDensityMatrix.h"
#include "NoPthreads.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "DiagBlockDiagMatrix.h"

namespace Dmrg {

template<typename TargetingType>
class DensityMatrixLocal : public DensityMatrixBase<TargetingType> {

	typedef DensityMatrixBase<TargetingType> BaseType;
	typedef typename TargetingType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename TargetingType::VectorWithOffsetType TargetVectorType;
	typedef typename TargetingType::TargetVectorType::value_type DensityMatrixElementType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef typename PsimagLite::Real<DensityMatrixElementType>::Type RealType;
	typedef typename DensityMatrixBase<TargetingType>::Params ParamsType;

public:

	typedef typename BaseType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BlockDiagonalMatrixType::BuildingBlockType BuildingBlockType;
	typedef ParallelDensityMatrix<BlockDiagonalMatrixType,
	BasisWithOperatorsType,
	TargetVectorType> ParallelDensityMatrixType;
	typedef PsimagLite::Parallelizer<ParallelDensityMatrixType> ParallelizerType;
	typedef typename TargetingType::VectorVectorVectorWithOffsetType
	VectorVectorVectorWithOffsetType;

	DensityMatrixLocal(const TargetingType& target,
	                   const LeftRightSuperType& lrs,
	                   const ParamsType& p)
	    :
	      progress_("DensityMatrixLocal"),
	      data_((p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.left() :
	                                                                            lrs.right()),
	      direction_(p.direction),
	      debug_(p.debug)
	{
		{
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Init partition for all targets";
			progress_.printline(msgg, std::cout);
		}

		const BasisWithOperatorsType& pBasis =
		        (p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.left() :
		                                                                        lrs.right();

		const BasisWithOperatorsType& pBasisSummed =
		        (p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.right() :
		                                                                        lrs.left();

		//loop over all partitions:
		for (SizeType m=0;m<pBasis.partition()-1;m++) {
			// size of this partition
			SizeType bs = pBasis.partition(m+1)-pBasis.partition(m);

			// density matrix block for this partition:
			BuildingBlockType matrixBlock(bs,bs);

			// weight of the ground state:
			RealType w = target.gsWeight();

			// if we are to target the ground state do it now:
			if (target.includeGroundStage()) {
				const VectorVectorVectorWithOffsetType& psi = target.psiConst();
				const SizeType nsectors = psi.size();

				for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
					const SizeType nexcited = psi[sectorIndex].size();

					for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {

						initPartition(matrixBlock,
						              pBasis,
						              m,
						              *(psi[sectorIndex][excitedIndex]),
						              pBasisSummed,
						              lrs.super(),
						              p.direction,
						              w);
					}
				}
			}

			// target all other states if any:
			for (SizeType ix = 0; ix < target.size(); ++ix) {
				RealType wnorm = target.normSquared(ix);
				if (fabs(wnorm) < 1e-6) continue;
				RealType w = target.weight(ix)/wnorm;
				initPartition(matrixBlock,
				              pBasis,
				              m,
				              target(ix),
				              pBasisSummed,
				              lrs.super(),
				              p.direction,
				              w);
			}

			// set this matrix block into data_
			data_.setBlock(m,pBasis.partition(m),matrixBlock);
		}

		{
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Done with init partition";
			progress_.printline(msgg, std::cout);
		}

	}

	virtual const BlockDiagonalMatrixType& operator()()
	{
		return data_;
	}

	void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz)
	{
		DiagBlockDiagMatrix<BlockDiagonalMatrixType>::diagonalise(data_,eigs,jobz);
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const DensityMatrixLocal& dm)
	{
		for (SizeType m = 0; m < dm.data_.blocks(); ++m) {
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<" ne="<<ne<<"\n";
			os<<dm.data_(m)<<"\n";
		}

		return os;
	}

private:

	void initPartition(BuildingBlockType& matrixBlock,
	                   BasisWithOperatorsType const &pBasis,
	                   SizeType m,
	                   const TargetVectorType& v,
	                   BasisWithOperatorsType const &pBasisSummed,
	                   BasisType const &pSE,
	                   ProgramGlobals::DirectionEnum direction,
	                   RealType weight)
	{
		ParallelDensityMatrixType helperDm(v,
		                                   pBasis,
		                                   pBasisSummed,
		                                   pSE,
		                                   direction,
		                                   m,
		                                   weight,
		                                   matrixBlock);
		ParallelizerType threadedDm(ConcurrencyType::codeSectionParams);
		threadedDm.loopCreate(helperDm);

	}

	ProgressIndicatorType progress_;
	BlockDiagonalMatrixType data_;
	ProgramGlobals::DirectionEnum direction_;
	bool debug_;
}; // class DensityMatrixLocal

} // namespace Dmrg

#endif

