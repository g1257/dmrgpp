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

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

	class ParallelPsiSplit {

	public:

		ParallelPsiSplit(const LeftRightSuperType& lrs,
		                 const VectorGenIjPatchType& vectorOfijPatches,
		                 const VectorWithOffset<ComplexOrRealType>& v,
		                 MatrixVectorType& allTargets)
		    : left_(lrs.left()),
		      right_(lrs.right()),
		      permInverse_(lrs.super().permutationInverse()),
		      v_(v),
		      allTargets_(allTargets)
		{
			ijPatch_ = vectorOfijPatches[vectorOfijPatches.size() - 1];
			for (SizeType i = 0; i < tasks(); ++i) {
				allTargets.push_back(new MatrixType());
			}
		}

		void doTask(SizeType ipatch, SizeType)
		{
			SizeType nl = left_.size();

			SizeType offset = v_.offset(0);
			SizeType igroup = ijPatch_->operator ()(GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = ijPatch_->operator ()(GenIjPatchType::RIGHT)[ipatch];

			SizeType sizeLeft = left_.partition(igroup+1) - left_.partition(igroup);
			SizeType sizeRight = right_.partition(jgroup+1) - right_.partition(jgroup);

			SizeType left_offset = left_.partition(igroup);
			SizeType right_offset = right_.partition(jgroup);

			MatrixType& matrix = *(allTargets_[ipatch]);
			matrix.resize(sizeLeft, sizeRight);

			for (SizeType ileft=0; ileft < sizeLeft; ++ileft) {
				for (SizeType iright=0; iright < sizeRight; ++iright) {

					SizeType i = ileft + left_offset;
					SizeType j = iright + right_offset;

					SizeType ij = i + j * nl;

					assert(i < nl);
					assert(j < right_.size());

					assert(ij < permInverse_.size());

					SizeType r = permInverse_[ij];
					if (r < offset || r >= offset + v_.effectiveSize(0))
						continue;

					matrix(ileft,iright) +=  v_.slowAccess(r);
				}
			}

		}

		SizeType tasks() const
		{
			return ijPatch_->operator ()(GenIjPatchType::LEFT).size();
		}

	private:

		const BasisWithOperatorsType& left_;
		const BasisWithOperatorsType& right_;
		const VectorSizeType& permInverse_;
		GenIjPatchType* ijPatch_;
		const VectorWithOffset<ComplexOrRealType>& v_;
		MatrixVectorType& allTargets_;
	};

	class ParallelSvd {

	public:

		ParallelSvd(const MatrixVectorType& allTargets,
		            MatrixType& mAll,
		            VectorRealType& eigs,
		            const LeftRightSuperType& lrs,
		            SizeType direction,
		            const VectorGenIjPatchType& vectorOfijPatches)
		    : allTargets_(allTargets),
		      mAll_(mAll),
		      eigs_(eigs),
		      lrs_(lrs),
		      direction_(direction)
		{
			ijPatch_ = vectorOfijPatches[vectorOfijPatches.size() - 1];
			SizeType oneSide = expandSys() ? lrs.left().size() : lrs.right().size();
			eigs.resize(oneSide, 0.0);
			mAll.resize(oneSide, oneSide);
		}

		void doTask(SizeType ipatch, SizeType)
		{
			MatrixType& m = *(allTargets_[ipatch]);
			SizeType freeSize = m.rows();
			MatrixType vt;
			VectorRealType eigsOnePatch(freeSize);

			svd('A', m, eigsOnePatch, vt);
			MatrixType* vMatrix = 0;
			if (!expandSys()) {
				vMatrix = new MatrixType();
				transposeConjugate(*vMatrix, vt);
			}

			saveThisPatch(mAll_, eigs_, m, vMatrix, eigsOnePatch, ipatch);
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
		                   SizeType ipatch)
		{
			const MatrixType& mLeftOrRight = expandSys() ? m : *vMatrix;
			const BasisType& basis = expandSys() ? lrs_.left() : lrs_.right();
			typename GenIjPatchType::LeftOrRightEnumType lOrR = (expandSys()) ?
			            GenIjPatchType::LEFT : GenIjPatchType::RIGHT;
			SizeType igroup = ijPatch_->operator ()(lOrR)[ipatch];
			SizeType offset = basis.partition(igroup);
			SizeType x = mLeftOrRight.rows();
			assert(x == mLeftOrRight.cols());

			for (SizeType i = 0; i < x; ++i) {
				eigs[i+offset] = eigsOnePatch[i];
				for (SizeType j = 0; j < x; ++j)
					mAll(i + offset, j + offset) += mLeftOrRight(i, j);
			}
		}

		bool expandSys() const
		{
			return (direction_ == ProgramGlobals::EXPAND_SYSTEM);
		}

		const MatrixVectorType& allTargets_;
		MatrixType& mAll_;
		VectorRealType& eigs_;
		const LeftRightSuperType& lrs_;
		SizeType direction_;
		GenIjPatchType* ijPatch_;
	};

public:

	DensityMatrixSvd(const TargettingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      params_(p),
	      lrs_(lrs),
	      vectorOfijPatches_(0)
	{
		{
			PsimagLite::OstringStream msg;
			msg<<"Init partition for all targets";
			progress_.printline(msg,std::cout);
		}

		SizeType oneOrZero = (target.includeGroundStage()) ? 1 : 0;
		SizeType targets = oneOrZero + target.size(); // Number of targets;

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
		for (SizeType i = 0; i < allTargets_.size(); ++i) {
			delete allTargets_[i];
			allTargets_[i] = 0;
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
		                        vectorOfijPatches_);
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
	                    const VectorWithOffset<ComplexOrRealType>& v,
	                    SizeType targets)
	{
		const BasisType& super = lrs_.super();

		SizeType m = v.sector(0);
		int state = super.partition(m);
		SizeType qn = super.qn(state);

		vectorOfijPatches_.push_back(new GenIjPatchType(lrs_, qn));
		typedef PsimagLite::Parallelizer<ParallelPsiSplit> ParallelizerType;
		ParallelizerType threaded(PsimagLite::Concurrency::npthreads,
		                          PsimagLite::MPI::COMM_WORLD);
		ParallelPsiSplit parallelPsiSplit(lrs_, vectorOfijPatches_, v, allTargets_);
		threaded.loopCreate(parallelPsiSplit);
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffsets<ComplexOrRealType>& v,
	                    SizeType targets)
	{
		err("useSvd doesn't yet work with VectorWithOffsets (sorry)\n");
	}

	ProgressIndicatorType progress_;
	const ParamsType& params_;
	MatrixVectorType allTargets_;
	SparseMatrixType data_;
	const LeftRightSuperType& lrs_;
	VectorGenIjPatchType vectorOfijPatches_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

