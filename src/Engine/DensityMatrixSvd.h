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
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
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
		                 const RealType& weight,
		                 SizeType targetsOffset)
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
		      weight_(weight),
		      targetsOffset_(targetsOffset)
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
			MatrixType& matrix = *(allTargets_[multiIndex + targetsOffset_*targetNumber_]);
			matrix.resize(sizeLeft, sizeRight);

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

					matrix(ileft,iright) +=  v_.slowAccess(r)*weight_;
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
		SizeType targetsOffset_;
	};

	class ParallelSvd {

		struct Opaque {
			Opaque(SizeType o, SizeType x1, SizeType r, SizeType c)
			    : offset(o), x(x1), rows(r), cols(c)
			{}

			SizeType offset;
			SizeType x;
			SizeType rows;
			SizeType cols;
		};

		typedef typename PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;
		typedef typename PsimagLite::Vector<Opaque*>::Type VectorOpaqueType;

	public:

		ParallelSvd(const MatrixVectorType& allTargets,
		            const LeftRightSuperType& lrs,
		            SizeType direction,
		            const VectorGenIjPatchType& vectorOfijPatches,
		            const VectorSizeType& patchBoundary,
		            SizeType targets,
		            SizeType targetsOffset)
		    : allTargets_(allTargets),
		      lrs_(lrs),
		      direction_(direction),
		      ijPatch_(vectorOfijPatches),
		      patchBoundary_(patchBoundary),
		      targetsOffset_(targetsOffset),
		      mAll_(targets, 0),
		      eigs_(targets, expandSys() ? lrs.left().size() : lrs.right().size()),
		      opaque_(patchBoundary_[patchBoundary_.size()-1], 0)
		{
			eigs_.setTo(0.0);
			SizeType oneSide = eigs_.cols();
			for (SizeType x = 0; x < mAll_.size(); ++x) {
				mAll_[x] = new MatrixType(oneSide, oneSide);
				mAll_[x]->setTo(0.0);
			}
		}

		~ParallelSvd()
		{
			for (SizeType x = 0; x < mAll_.size(); ++x) {
				delete mAll_[x];
				mAll_[x] = 0;
			}

			for (SizeType x = 0; x < opaque_.size(); ++x) {
				delete opaque_[x];
				opaque_[x] = 0;
			}
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

			saveThisPatch(m, vMatrix, eigsOnePatch, multiIndex);
			delete vMatrix;
		}

		SizeType tasks() const
		{
			return allTargets_.size();
		}

		void sync(SparseMatrixType& data, VectorRealType& eigs) const
		{
			eigs.resize(eigs_.cols(), 0.0);
			MatrixType m(mAll_[0]->rows(), mAll_[0]->cols());
			m.setTo(0.0);
			for (SizeType i = 0; i < ijPatch_.size(); ++i) {
				const GenIjPatchType& ijPatch = *(ijPatch_[i]);
				for (SizeType j = 0; j < ijPatch(GenIjPatchType::LEFT).size(); ++j) {
					SizeType reducedIndex = patchBoundary_[i] + j;
					SizeType t = findLowestEigIndex(reducedIndex);
					copyOne(m, eigs, reducedIndex, t);
				}
			}

			fullMatrixToCrsMatrix(data, m);
		}

	private:

		void saveThisPatch(const MatrixType& m,
		                   const MatrixType* vMatrix,
		                   const VectorRealType& eigsOnePatch,
		                   SizeType multiIndex)
		{
			const MatrixType& mLeftOrRight = expandSys() ? m : *vMatrix;
			const BasisType& basis = expandSys() ? lrs_.left() : lrs_.right();
			typename GenIjPatchType::LeftOrRightEnumType lOrR = (expandSys()) ?
			            GenIjPatchType::LEFT : GenIjPatchType::RIGHT;
			SizeType t = static_cast<SizeType>(multiIndex/targetsOffset_);
			assert(t < mAll_.size());
			SizeType reducedIndex = multiIndex % targetsOffset_;
			PairSizeType p = multiIndexToPatch(reducedIndex);
			SizeType igroup = ijPatch_[p.first]->operator()(lOrR)[p.second];
			SizeType offset = basis.partition(igroup);
			SizeType partSize = basis.partition(igroup + 1) - offset;
			SizeType rows = mLeftOrRight.rows();
			SizeType cols = mLeftOrRight.cols();

			for (SizeType i = 0; i < rows; ++i)
				for (SizeType j = 0; j < cols; ++j)
					mAll_[t]->operator()(i + offset, j + offset) += mLeftOrRight(i, j);

			SizeType x = eigsOnePatch.size();
			if (x > partSize) x = partSize;
			for (SizeType i = 0; i < x; ++i)
				eigs_(t, i + offset) = eigsOnePatch[i]*eigsOnePatch[i];

			assert(reducedIndex < opaque_.size());
			if (opaque_[reducedIndex] != 0) return;
			opaque_[reducedIndex] = new Opaque(offset,  x, rows, cols);
		}

		void copyOne(MatrixType& mAll,
		             VectorRealType& eigs,
		             SizeType reducedIndex,
		             SizeType t) const
		{
			assert(t < mAll_.size());
			SizeType offset = opaque_[reducedIndex]->offset;
			SizeType rows = opaque_[reducedIndex]->rows;
			SizeType cols = opaque_[reducedIndex]->cols;

			for (SizeType i = 0; i < rows; ++i)
				for (SizeType j = 0; j < cols; ++j)
					mAll(i + offset, j + offset) = mAll_[t]->operator()(i + offset, j + offset);

			SizeType x = opaque_[reducedIndex]->x;
			assert(x + offset <= eigs.size());
			for (SizeType i = 0; i < x; ++i)
				eigs[i + offset] = eigs_(t, i + offset);
		}

		SizeType findLowestEigIndex(SizeType reducedIndex) const
		{
			// FIXME: Precompute sum of eigs for each target
			SizeType tmin = 0;
			SizeType min = sumOfEigs(tmin, reducedIndex);
			for (SizeType t = 1; t < mAll_.size(); ++t) {
				RealType sum = sumOfEigs(t, reducedIndex);
				if (sum < min) {
					min = sum;
					tmin = t;
				}
			}

			return tmin;
		}

		RealType sumOfEigs(SizeType t, SizeType reducedIndex) const
		{
			SizeType x = opaque_[reducedIndex]->x;
			SizeType offset = opaque_[reducedIndex]->offset;

			RealType sum = 0.0;
			for (SizeType i = 0; i < x; ++i)
				sum += eigs_(t, i + offset);

			return sum;
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
		const LeftRightSuperType& lrs_;
		SizeType direction_;
		const VectorGenIjPatchType& ijPatch_;
		const VectorSizeType& patchBoundary_;
		SizeType targetsOffset_;
		VectorMatrixType mAll_;
		MatrixRealType eigs_;
		VectorOpaqueType opaque_;
	};

public:

	DensityMatrixSvd(const TargettingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      params_(p),
	      lrs_(lrs),
	      targets_(target.size() + (target.includeGroundStage()) ? 1 : 0),
	      targetOffset_(0)
	{
		for (SizeType x = 0; x  < targets_; ++x) {
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
				patchBoundary_.push_back(targetOffset_);
				qnToPatch_.push_back(qn);
				targetOffset_ += ptr->operator ()(GenIjPatchType::LEFT).size();
			}
		}

		patchBoundary_.push_back(targetOffset_);
		allTargets_.resize(targetOffset_*targets_);
		for (SizeType i = 0; i < allTargets_.size(); ++i)
			allTargets_[i] = new MatrixType();

		{
			PsimagLite::OstringStream msg;
			msg<<"Found "<<allTargets_.size()<<" multindices, targetOffset= ";
			msg<<targetOffset_<<", qns= "<<qnToPatch_.size();
			progress_.printline(msg,std::cout);
		}

		RealType weights = 0.0;
		for (SizeType x = 0; x < targets_; ++x)
			weights += addThisTarget(x, target);

		{
			PsimagLite::OstringStream msg;
			msg<<"Done with init partition, targets= "<<targets_;
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
		typedef PsimagLite::Parallelizer<ParallelSvd> ParallelizerType;
		ParallelizerType threaded(PsimagLite::Concurrency::npthreads,
		                          PsimagLite::MPI::COMM_WORLD);
		ParallelSvd parallelSvd(allTargets_,
		                        lrs_,
		                        params_.direction,
		                        vectorOfijPatches_,
		                        patchBoundary_,
		                        targets_,
		                        targetOffset_);
		threaded.loopCreate(parallelSvd);

		parallelSvd.sync(data_, eigs);
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
	                       const TargettingType& target)

	{
		SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

		const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
		            target.gs() : target(x2);

		RealType weight = (target.includeGroundStage() && x == 0 ) ?
		            target.gsWeight() : target.weight(x2);

		addThisTarget2(x, v, sqrt(weight));
		return weight;
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffsetType& v,
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
			                                  targets_,
			                                  weight,
			                                  targetOffset_);
			threaded.loopCreate(parallelPsiSplit);
		}
	}

	ProgressIndicatorType progress_;
	const ParamsType& params_;
	MatrixVectorType allTargets_;
	SparseMatrixType data_;
	const LeftRightSuperType& lrs_;
	SizeType targets_;
	SizeType targetOffset_;
	VectorGenIjPatchType vectorOfijPatches_;
	VectorSizeType qnToPatch_;
	VectorSizeType patchBoundary_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

