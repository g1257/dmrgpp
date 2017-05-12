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

#ifndef DENSITY_MATRIX_SVD_H
#define DENSITY_MATRIX_SVD_H
#include "ProgressIndicator.h"
#include "TypeToString.h"
#include "BlockMatrix.h"
#include "DensityMatrixBase.h"
#include "NoPthreads.h"
#include "Concurrency.h"
#include "Parallelizer.h"

namespace Dmrg {

template<typename TargettingType>
class DensityMatrixSvd : public DensityMatrixBase<TargettingType> {

	typedef DensityMatrixBase<TargettingType> BaseType;
	typedef typename TargettingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename BaseType::BuildingBlockType MatrixType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename BaseType::Params ParamsType;

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

public:

	typedef typename BaseType::BlockMatrixType BlockMatrixType;
	typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;
	typedef ParallelDensityMatrix<BlockMatrixType,
	BasisWithOperatorsType,
	VectorWithOffsetType> ParallelDensityMatrixType;
	typedef PsimagLite::Parallelizer<ParallelDensityMatrixType> ParallelizerType;

	DensityMatrixSvd(const TargettingType&,
	                   const BasisWithOperatorsType& pBasis,
	                   const BasisWithOperatorsType&,
	                   const BasisType&,
	                   const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      debug_(p.debug),
	      verbose_(p.verbose)
	{}

	virtual SparseMatrixType& operator()()
	{
		return data_;
	}

	virtual SizeType rank() { return allTargets_.rows(); }

	void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz)
	{
		err("Svd not implemented\n");
		fullMatrixToCrsMatrix(data_, allTargets_);
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

		SizeType oneOrZero = (target.includeGroundStage()) ? 1 : 0;
		SizeType targets = oneOrZero + target.size(); // Number of targets;
		SizeType freeSize = pBasis.size();
		SizeType summedSize = pBasisSummed.size();
		allTargets_.resize(targets*summedSize, freeSize);
		for (SizeType x = 0; x < targets; ++x)
			addThisTarget(x, freeSize, summedSize, target, p.direction, pSE, targets);

		{
			PsimagLite::OstringStream msg;
			msg<<"Done with init partition";
			progress_.printline(msg,std::cout);
		}
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const DensityMatrixSvd& dm)
	{
		for (SizeType m = 0; m < dm.data_.blocks(); ++m) {
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<" ne="<<ne<<"\n";
			os<<dm.data_(m)<<"\n";
		}

		return os;
	}

private:

	void addThisTarget(SizeType x,
	                   SizeType freeSize,
	                   SizeType summedSize,
	                   const TargettingType& target,
	                   SizeType direction,
	                   const BasisType& pSE,
	                   SizeType targets)
	{
		SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

		const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
		            target.gs() : target(x2);

		for (SizeType alpha = 0; alpha < freeSize; ++alpha) {
			for (SizeType beta = 0; beta < summedSize; ++beta) {
				SizeType ind = (direction == ProgramGlobals::EXPAND_SYSTEM) ?
				            alpha + beta*summedSize : beta + alpha*freeSize;
				ind = pSE.permutationInverse(ind);
				allTargets_(x + beta*targets, alpha) = v.slowAccess(ind);
			}
		}
	}

	ProgressIndicatorType progress_;
	MatrixType allTargets_;
	SparseMatrixType data_;
	bool debug_;
	bool verbose_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

