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
#include "DensityMatrixBase.h"
#include "NoPthreads.h"
#include "Concurrency.h"
#include "MatrixVectorKron/GenIjPatch.h"

namespace Dmrg {

template<typename TargettingType>
class DensityMatrixSvd : public DensityMatrixBase<TargettingType> {

	typedef DensityMatrixBase<TargettingType> BaseType;
	typedef typename TargettingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename TargettingType::LeftRightSuperType LeftRightSuperType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename TargettingType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename BaseType::BuildingBlockType MatrixType;
	typedef PsimagLite::Matrix<SizeType> MatrixSizeType;
	typedef typename PsimagLite::Vector<MatrixType*>::Type MatrixVectorType;
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename BaseType::Params ParamsType;
	typedef GenIjPatch<LeftRightSuperType> GenIjPatchType;
	typedef typename GenIjPatchType::VectorSizeType VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<GenIjPatchType*>::Type VectorGenIjPatchType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef typename BaseType::BlockDiagonalMatrixType BlockDiagonalMatrixType;

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

	class GroupsStruct {

	};

	typedef GroupsStruct GroupsStructType;

	class ParallelPsiSplit {

	public:

		ParallelPsiSplit(const LeftRightSuperType& lrs,
		                 const GenIjPatchType& ijPatch,
		                 const VectorWithOffsetType& v,
		                 SizeType target,
		                 SizeType sector,
		                 MatrixVectorType& allTargets)
		    : lrs_(lrs),
		      ijPatch_(ijPatch),
		      v_(v),
		      target_(target),
		      sector_(sector),
		      allTargets_(allTargets)
		{}

		void doTask(SizeType ipatch, SizeType)
		{
			SizeType igroup = ijPatch_(GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = ijPatch_(GenIjPatchType::RIGHT)[ipatch];

			SizeType groupToTruncate = (allTargets_.direction() == ProgramGlobals::EXPAND_SYSTEM) ?
			            igroup : jgroup;

			packMatrix(groupToTruncate);
		}

		SizeType tasks() const
		{
			return ijPatch_(GenIjPatchType::LEFT).size();
		}

	private:

		void packMatrix(SizeType igroup)
		{
			assert(group < allTargets_.size());
			MatrixType& matrix = *(allTargets_[igroup]);
			SizeType m = v_.sector(sector_);
			SizeType offset = v_.offset(m);
			SizeType sizeLeft = lrs_.left().partition(igroup+1) - lrs_.left().partition(igroup);
			SizeType sizeRight = lrs_.right().partition(jgroup+1) - lrs_.right().partition(jgroup);
			matrix.resize(sizeLeft, sizeRight);

			SizeType nl = lrs_.left().size();
			SizeType leftOffset = lrs_.left().partition(igroup);
			SizeType rightOffset = lrs_.right().partition(jgroup);

			for (SizeType ileft=0; ileft < sizeLeft; ++ileft) {
				for (SizeType iright=0; iright < sizeRight; ++iright) {

					SizeType i = ileft + leftOffset;
					SizeType j = iright + rightOffset;

					SizeType ij = i + j * nl;

					assert(i < nl);
					assert(j < lrs_.right().size());

					assert(ij < lrs_.super().permutationInverse().size());

					SizeType r = lrs_.super().permutationInverse()[ij];
					if (r < offset || r >= offset + v_.effectiveSize(m))
						continue;

					matrix(ileft,iright + additionalOffset_) +=  v_.slowAccess(r);
				}
			}
		}

		SizeType ipatchToMultiIndex(SizeType qn, SizeType ipatch) const
		{
			typename VectorSizeType::const_iterator it = std::find(qnToPatch_.begin(),
			                                                       qnToPatch_.end(),
			                                                       qn);
			assert(it != qnToPatch_.end());
			SizeType ind = it - qnToPatch_.begin();
			return patchBoundary_[ind] + ipatch;
		}

		const LeftRightSuperType& lrs_;
		const GenIjPatchType& ijPatch_;
		const VectorWithOffsetType& v_;
		SizeType target_;
		SizeType sector_;
		GroupsStructType& allTargets_;
	};

	class ParallelSvd {

	public:

		ParallelSvd(const GroupsStructType& allTargets,
		            VectorRealType& eigs)
		    : allTargets_(allTargets),
		      eigs_(eigs)
		{}

		void doTask(SizeType ipatch, SizeType)
		{
			MatrixType& m = allTargets_.matrix(ipatch);
			MatrixType vt;
			VectorRealType eigsOnePatch;

			svd('S', m, eigsOnePatch, vt);

			const BasisType& basis = allTargets_.basis();

			SizeType igroup = allTargets_.groupIndex(ipatch);
			SizeType offset = basis.partition(igroup);
			SizeType partSize = basis.partition(igroup + 1) - offset;

			SizeType x = eigsOnePatch.size();
			if (x > partSize) x = partSize;
			assert(x + offset <= eigs.size());
			for (SizeType i = 0; i < x; ++i)
				eigs[i + offset] = eigsOnePatch[i]*eigsOnePatch[i];
		}

		SizeType tasks() const
		{
			return allTargets_.size();
		}

		const BlockDiagonalMatrixType& blockMatrix() const
		{
			return *mAll_;
		}

	private:

		const MatrixVectorType& allTargets_;
		VectorRealType& eigs_;
	};

public:

	DensityMatrixSvd(const TargettingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      lrs_(lrs),
	      params_(p),
	      data_(0),
	      allTargets_(lrs, p.direction)
	{
		SizeType oneOrZero = (target.includeGroundStage()) ? 1 : 0;
		SizeType targets = oneOrZero + target.size(); // Number of targets;

		for (SizeType x = 0; x  < targets; ++x) {
			SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

			const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
			            target.gs() : target(x2);

			SizeType sectors = v.sectors();
			for (SizeType sector = 0; sector < sectors; ++sector) {	
				GenIjPatchType* ptr = new GenIjPatchType(lrs, qn);
				vectorOfijPatches_.push_back(ptr);
				const VectorSizeType& groups =  ptr->operator()(direction);
				for (SizeType i = 0; i < groups.size(); ++i) {
					SizeType igroup = groups[i];
					SizeType jgroup = ptr->operator()(directionPrime)[igroup];
					// have I seen this group before? No --> create group
					if (it == seenGroups.end()) {
						seenGroups.push_back(igroup);
						allTargets_.push(igroup,jgroup);
					} else { //  Yes, add to group
						SizeType groupIndex = it - seenGroups.begin();
						assert(groupIndex < groupsStruct.size());
						allTargets_.addToGroup(groupIndex, jgroup);
					}
				}
			}
		}

		allTargets_.finalize();

		{
			PsimagLite::OstringStream msg;
			msg<<"Found "<<groupsIncluded.size()<<" groups on ";
			progress_.printline(msg,std::cout);
		}

		for (SizeType x = 0; x < targets; ++x)
			addThisTarget(x, target, targets);

		{
			PsimagLite::OstringStream msg;
			msg<<"Done with init partition";
			progress_.printline(msg,std::cout);
		}
	}

	~DensityMatrixSvd()
	{
		for (SizeType i = 0; i < vectorOfijPatches_.size(); ++i) {
			delete vectorOfijPatches_[i];
			vectorOfijPatches_[i] = 0;
		}

		delete data_;
		data_ = 0;
	}

	virtual const BlockDiagonalMatrixType& operator()()
	{
		return *data_;
	}

	void diag(VectorRealType& eigs,char jobz)
	{
		typedef PsimagLite::Parallelizer<ParallelSvd> ParallelizerType;
		ParallelizerType threaded(PsimagLite::Concurrency::npthreads,
		                          PsimagLite::MPI::COMM_WORLD);
		ParallelSvd parallelSvd(allTargets_,
		                        eigs);
		threaded.loopCreate(parallelSvd);

		data_ = &(parallelSvd.blockMatrix());
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const DensityMatrixSvd& dm)
	{
		for (SizeType m = 0; m < dm.data_->blocks(); ++m) {
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<" ne="<<ne<<"\n";
			os<<dm.data_->operator()(m)<<"\n";
		}

		return os;
	}

private:

	void addThisTarget(SizeType x,
	                   const TargettingType& target,
	                   SizeType targets)

	{
		SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

		const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
		            target.gs() : target(x2);

		addThisTarget2(x, v, targets);
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffsetType& v,
	                    SizeType targets)
	{
		const BasisType& super = lrs_.super();
		assert(targets == 1 && x == 0);
		for (SizeType sector = 0; sector < v.sectors(); ++sector) {
			SizeType m = v.sector(sector);
			int state = super.partition(m);
			SizeType qn = super.qn(state);

			typename VectorSizeType::iterator it = std::find(qnToPatch_.begin(),
			                                                 qnToPatch_.end(),
			                                                 qn);
			assert(it != qnToPatch_.end());
			SizeType ind = it - qnToPatch_.begin();
			typedef PsimagLite::Parallelizer<ParallelPsiSplit> ParallelizerType;
			ParallelizerType threaded(PsimagLite::Concurrency::npthreads,
			                          PsimagLite::MPI::COMM_WORLD);
			ParallelPsiSplit parallelPsiSplit(lrs_,
			                                  *(vectorOfijPatches_[ind]),
			                                  v,
			                                  x,
			                                  sector,
			                                  allTargets_);
			threaded.loopCreate(parallelPsiSplit);
		}
	}

	ProgressIndicatorType progress_;
	const LeftRightSuperType& lrs_;
	const ParamsType& params_;
	const BlockDiagonalMatrixType* data_;
	GroupsStructType allTargets_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

