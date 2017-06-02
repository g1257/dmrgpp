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

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

	class ParallelPsiSplit {

	public:

		ParallelPsiSplit(const LeftRightSuperType& lrs,
		                 const GenIjPatchType& ijPatch,
		                 const VectorWithOffsetType& v,
		                 SizeType sector,
		                 MatrixVectorType& allTargets,
		                 const VectorSizeType& qnToPatch,
		                 const VectorSizeType& patchBoundary,
		                 SizeType direction,
		                 SizeType targetNumber,
		                 SizeType targets,
		                 const RealType& weight)
		    : lrs_(lrs),
		      ijPatch_(ijPatch),
		      v_(v),
		      sector_(sector),
		      allTargets_(allTargets),
		      qnToPatch_(qnToPatch),
		      patchBoundary_(patchBoundary),
		      direction_(direction),
		      targetNumber_(targetNumber),
		      targets_(targets),
		      weight_(weight)
		{}

		void doTask(SizeType ipatch, SizeType)
		{
			SizeType nl = lrs_.left().size();

			SizeType m = v_.sector(sector_);
			SizeType offset = v_.offset(m);
			SizeType igroup = ijPatch_(GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = ijPatch_(GenIjPatchType::RIGHT)[ipatch];

			SizeType sizeLeft = lrs_.left().partition(igroup+1) - lrs_.left().partition(igroup);
			SizeType sizeRight = lrs_.right().partition(jgroup+1) - lrs_.right().partition(jgroup);
			SizeType leftOffset = lrs_.left().partition(igroup);
			SizeType rightOffset = lrs_.right().partition(jgroup);

			int state = lrs_.super().partition(m);
			SizeType qn = lrs_.super().qn(state);
			SizeType multiIndex = ipatchToMultiIndex(qn, ipatch);
			MatrixType& matrix = *(allTargets_[multiIndex]);
			matrix.resize(expandSys() ? sizeLeft : sizeLeft*targets_,
			              expandSys() ? sizeRight*targets_ : sizeRight);

			for (SizeType ileft=0; ileft < sizeLeft; ++ileft) {
				SizeType ileftModif = (expandSys()) ?
				            ileft : ileft + targetNumber_*sizeLeft;
				for (SizeType iright=0; iright < sizeRight; ++iright) {
					SizeType irightModif = (expandSys()) ?
					            iright + targetNumber_*sizeRight : iright;
					SizeType i = ileft + leftOffset;
					SizeType j = iright + rightOffset;

					SizeType ij = i + j * nl;

					assert(i < nl);
					assert(j < lrs_.right().size());

					assert(ij < lrs_.super().permutationInverse().size());

					SizeType r = lrs_.super().permutationInverse()[ij];
					if (r < offset || r >= offset + v_.effectiveSize(m))
						continue;

					matrix(ileftModif,irightModif) +=  v_.slowAccess(r)*weight_;
				}
			}
		}

		SizeType tasks() const
		{
			return ijPatch_(GenIjPatchType::LEFT).size();
		}

	private:

		SizeType ipatchToMultiIndex(SizeType qn, SizeType ipatch) const
		{
			typename VectorSizeType::const_iterator it = std::find(qnToPatch_.begin(),
			                                                       qnToPatch_.end(),
			                                                       qn);
			assert(it != qnToPatch_.end());
			SizeType ind = it - qnToPatch_.begin();
			return patchBoundary_[ind] + ipatch;
		}

		bool expandSys() const
		{
			return (direction_ == ProgramGlobals::EXPAND_SYSTEM);
		}

		const LeftRightSuperType& lrs_;
		const GenIjPatchType& ijPatch_;
		const VectorWithOffsetType& v_;
		SizeType sector_;
		MatrixVectorType& allTargets_;
		const VectorSizeType& qnToPatch_;
		const VectorSizeType& patchBoundary_;
		SizeType direction_;
		SizeType targetNumber_;
		SizeType targets_;
		RealType weight_;
	};

	class ParallelSvd {

	public:

		ParallelSvd(const MatrixVectorType& allTargets,
		            MatrixType& mAll,
		            VectorRealType& eigs,
		            const LeftRightSuperType& lrs,
		            SizeType direction,
		            const VectorGenIjPatchType& vectorOfijPatches,
		            const VectorSizeType& patchBoundary)
		    : allTargets_(allTargets),
		      mAll_(mAll),
		      eigs_(eigs),
		      lrs_(lrs),
		      direction_(direction),
		      ijPatch_(vectorOfijPatches),
		      patchBoundary_(patchBoundary)
		{
			SizeType oneSide = expandSys() ? lrs.left().size() : lrs.right().size();
			eigs_.resize(oneSide);
			std::fill(eigs_.begin(), eigs_.end(), 0.0);
			mAll.resize(oneSide, oneSide);
		}

		void doTask(SizeType multiIndex, SizeType)
		{
			MatrixType& m = *(allTargets_[multiIndex]);
			MatrixType vt;
			VectorRealType eigsOnePatch;

			svd('S', m, eigsOnePatch, vt);
			MatrixType* vMatrix = 0;
			if (!expandSys()) {
				vMatrix = new MatrixType();
				transposeConjugate(*vMatrix, vt);
			}

			saveThisPatch(mAll_, eigs_, m, vMatrix, eigsOnePatch, multiIndex);
			delete vMatrix;
		}

		SizeType tasks() const
		{
			return allTargets_.size();
		}

	private:

		void saveThisPatch(MatrixType& mAll,
		                   VectorRealType& eigs,
		                   const MatrixType& m,
		                   const MatrixType* vMatrix,
		                   const VectorRealType& eigsOnePatch,
		                   SizeType multiIndex)
		{
			const MatrixType& mLeftOrRight = expandSys() ? m : *vMatrix;
			const BasisType& basis = expandSys() ? lrs_.left() : lrs_.right();
			typename GenIjPatchType::LeftOrRightEnumType lOrR = (expandSys()) ?
			            GenIjPatchType::LEFT : GenIjPatchType::RIGHT;
			PairSizeType p = multiIndexToPatch(multiIndex);
			SizeType igroup = ijPatch_[p.first]->operator ()(lOrR)[p.second];
			SizeType offset = basis.partition(igroup);
			SizeType partSize = basis.partition(igroup + 1) - offset;
			SizeType rows = mLeftOrRight.rows();
			SizeType cols = mLeftOrRight.cols();

			for (SizeType i = 0; i < rows; ++i)
				for (SizeType j = 0; j < cols; ++j)
					mAll(i + offset, j + offset) += mLeftOrRight(i, j);

			SizeType x = eigsOnePatch.size();
			if (x > partSize) x = partSize;
			assert(x + offset <= eigs.size());
			for (SizeType i = 0; i < x; ++i)
				eigs[i + offset] = eigsOnePatch[i]*eigsOnePatch[i];

		}

		bool expandSys() const
		{
			return (direction_ == ProgramGlobals::EXPAND_SYSTEM);
		}

		PairSizeType multiIndexToPatch(SizeType multiIndex) const
		{
			assert(patchBoundary_.size() > 1);
			for (SizeType i = 0; i < patchBoundary_.size(); ++i) {
				if (multiIndex >= patchBoundary_[i] && multiIndex < patchBoundary_[i+1])
					return PairSizeType(i, multiIndex - patchBoundary_[i]);
			}

			err("DensityMatrixSvd: multiIndexToPatch internal error\n");
			throw PsimagLite::RuntimeError("unreachable\n");
		}

		const MatrixVectorType& allTargets_;
		MatrixType& mAll_;
		VectorRealType& eigs_;
		const LeftRightSuperType& lrs_;
		SizeType direction_;
		const VectorGenIjPatchType& ijPatch_;
		const VectorSizeType& patchBoundary_;
	};

public:

	DensityMatrixSvd(const TargettingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      params_(p),
	      lrs_(lrs)
	{
		SizeType oneOrZero = (target.includeGroundStage()) ? 1 : 0;
		SizeType targets = oneOrZero + target.size(); // Number of targets;

		SizeType sum = 0;
		for (SizeType x = 0; x  < targets; ++x) {
			SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

			const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
			            target.gs() : target(x2);

			SizeType sectors = v.sectors();
			for (SizeType sector = 0; sector < sectors; ++sector) {
				SizeType m = v.sector(sector);
				int state = lrs.super().partition(m);
				SizeType qn = lrs.super().qn(state);
				if (std::find(qnToPatch_.begin(), qnToPatch_.end(), qn) != qnToPatch_.end())
					continue;
				GenIjPatchType* ptr = new GenIjPatchType(lrs, qn);
				vectorOfijPatches_.push_back(ptr);
				patchBoundary_.push_back(sum);
				qnToPatch_.push_back(qn);
				sum += ptr->operator ()(GenIjPatchType::LEFT).size();
			}
		}

		patchBoundary_.push_back(sum);
		allTargets_.resize(sum);
		for (SizeType i = 0; i < sum; ++i)
			allTargets_[i] = new MatrixType();

		{
			PsimagLite::OstringStream msg;
			msg<<"Found "<<sum<<" multindices and "<<qnToPatch_.size()<<" qns";
			progress_.printline(msg,std::cout);
		}

		RealType weights = 0.0;
		for (SizeType x = 0; x < targets; ++x)
			weights += addThisTarget(x, target, targets);

		{
			PsimagLite::OstringStream msg;
			msg<<"Done with init partition, targets= "<<targets;
			msg<<" sum of weights= "<<weights;
			progress_.printline(msg,std::cout);
		}
	}

	~DensityMatrixSvd()
	{
		for (SizeType i = 0; i < allTargets_.size(); ++i) {
			delete allTargets_[i];
			allTargets_[i] = 0;
		}

		for (SizeType i = 0; i < vectorOfijPatches_.size(); ++i) {
			delete vectorOfijPatches_[i];
			vectorOfijPatches_[i] = 0;
		}
	}

	virtual SparseMatrixType& operator()()
	{
		return data_;
	}

	void diag(VectorRealType& eigs,char jobz)
	{
		MatrixType mAll;
		typedef PsimagLite::Parallelizer<ParallelSvd> ParallelizerType;
		ParallelizerType threaded(PsimagLite::Concurrency::npthreads,
		                          PsimagLite::MPI::COMM_WORLD);
		ParallelSvd parallelSvd(allTargets_,
		                        mAll,
		                        eigs,
		                        lrs_,
		                        params_.direction,
		                        vectorOfijPatches_,
		                        patchBoundary_);
		threaded.loopCreate(parallelSvd);

		fullMatrixToCrsMatrix(data_, mAll);
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

	RealType addThisTarget(SizeType x,
	                       const TargettingType& target,
	                       SizeType targets)

	{
		SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

		const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
		            target.gs() : target(x2);

		RealType weight = (target.includeGroundStage() && x == 0 ) ?
		            target.gsWeight() : target.weight(x2);

		addThisTarget2(x, v, targets, sqrt(weight));
		return weight;
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffsetType& v,
	                    SizeType targets,
	                    const RealType& weight)
	{
		const BasisType& super = lrs_.super();
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
			                                  sector,
			                                  allTargets_,
			                                  qnToPatch_,
			                                  patchBoundary_,
			                                  params_.direction,
			                                  x,
			                                  targets,
			                                  weight);
			threaded.loopCreate(parallelPsiSplit);
		}
	}

	ProgressIndicatorType progress_;
	const ParamsType& params_;
	MatrixVectorType allTargets_;
	SparseMatrixType data_;
	const LeftRightSuperType& lrs_;
	VectorGenIjPatchType vectorOfijPatches_;
	VectorSizeType qnToPatch_;
	VectorSizeType patchBoundary_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

