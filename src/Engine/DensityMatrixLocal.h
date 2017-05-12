
/*
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
#include "BlockMatrix.h"
#include "DensityMatrixBase.h"
#include "ParallelDensityMatrix.h"
#include "NoPthreads.h"
#include "Concurrency.h"
#include "Parallelizer.h"

namespace Dmrg {

template<typename TargettingType>
class DensityMatrixLocal : public DensityMatrixBase<TargettingType> {

	typedef DensityMatrixBase<TargettingType> BaseType;
	typedef typename TargettingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename TargettingType::VectorWithOffsetType TargetVectorType;
	typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef typename PsimagLite::Real<DensityMatrixElementType>::Type RealType;
	typedef typename DensityMatrixBase<TargettingType>::Params ParamsType;

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

public:

	typedef typename BaseType::BlockMatrixType BlockMatrixType;
	typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;
	typedef ParallelDensityMatrix<BlockMatrixType,
	BasisWithOperatorsType,
	TargetVectorType> ParallelDensityMatrixType;
	typedef PsimagLite::Parallelizer<ParallelDensityMatrixType> ParallelizerType;

	DensityMatrixLocal(const TargettingType&,
	                   const BasisWithOperatorsType& pBasis,
	                   const BasisWithOperatorsType&,
	                   const BasisType&,
	                   const ParamsType& p)
	    :
	      progress_("DensityMatrixLocal"),
	      data_(pBasis.size(),
	            pBasis.partition()-1),
	      direction_(p.direction),
	      debug_(p.debug),
	      verbose_(p.verbose)
	{}

	virtual BlockMatrixType& operator()()
	{
		return data_;
	}

	virtual SizeType rank() { return data_.rank(); }

	void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz)
	{
		diagonalise(data_,eigs,jobz);
	}

	virtual void init(const TargettingType& target,
	                  BasisWithOperatorsType const &pBasis,
	                  const BasisWithOperatorsType& pBasisSummed,
	                  BasisType const &pSE,
	                  const ParamsType& p)
	{
		{
			PsimagLite::OstringStream msg;
			msg<<"Init partition for all targets";
			progress_.printline(msg,std::cout);
		}

		//loop over all partitions:
		for (SizeType m=0;m<pBasis.partition()-1;m++) {
			// size of this partition
			SizeType bs = pBasis.partition(m+1)-pBasis.partition(m);

			// density matrix block for this partition:
			BuildingBlockType matrixBlock(bs,bs);

			// weight of the ground state:
			RealType w = target.gsWeight();

			// if we are to target the ground state do it now:
			if (target.includeGroundStage()) initPartition(matrixBlock,
			                                               pBasis,
			                                               m,
			                                               target.gs(),
			                                               pBasisSummed,
			                                               pSE,
			                                               p.direction,
			                                               w);

			// target all other states if any:
			for (SizeType ix = 0; ix < target.size(); ++ix) {
				RealType wnorm = target.normSquared(ix);
				if (fabs(wnorm) < 1e-6) continue;
				RealType w = target.weight(ix)/wnorm;
				initPartition(matrixBlock,pBasis,m,target(ix),
				              pBasisSummed,pSE,p.direction,w);
			}

			// set this matrix block into data_
			data_.setBlock(m,pBasis.partition(m),matrixBlock);
		}
		{
			PsimagLite::OstringStream msg;
			msg<<"Done with init partition";
			progress_.printline(msg,std::cout);
		}
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
	                   SizeType direction,
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
		ParallelizerType threadedDm(ConcurrencyType::npthreads,
		                            PsimagLite::MPI::COMM_WORLD);
		threadedDm.loopCreate(helperDm);

	}

	ProgressIndicatorType progress_;
	BlockMatrixType data_;
	SizeType direction_;
	bool debug_;
	bool verbose_;
}; // class DensityMatrixLocal

} // namespace Dmrg

#endif

