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
	typedef typename GenIjPatchType::GenGroupType GenGroupType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<GenIjPatchType*>::Type VectorGenIjPatchType;

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

public:

	DensityMatrixSvd(const TargettingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      params_(p),
	      lrs_(lrs),
	      gengroupLeft_(lrs.left()),
	      gengroupRight_(lrs.right()),
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
		SizeType npatches = allTargets_.size();
		SizeType oneSide = expandSys() ? lrs_.left().size() : lrs_.right().size();
		eigs.resize(oneSide, 0.0);
		MatrixType mAll(oneSide, oneSide);
		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {
			MatrixType& m = *(allTargets_[ipatch]);
			SizeType freeSize = m.rows();
			MatrixType vt(freeSize, freeSize);
			VectorRealType eigsOnePatch(freeSize);

			svd('A', m, eigsOnePatch, vt);
			saveThisPatch(mAll, eigs, m, vt, eigsOnePatch, ipatch);
		}

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
		const BasisWithOperatorsType& left = lrs_.left();
		const BasisWithOperatorsType& right = lrs_.right();

		SizeType m = v.sector(0);
		int state = super.partition(m);
		SizeType qn = super.qn(state);

		vectorOfijPatches_.push_back(new GenIjPatchType(lrs_, qn));
		GenIjPatchType& ijPatch = *(vectorOfijPatches_[vectorOfijPatches_.size() - 1]);

		const VectorSizeType& permInverse = super.permutationInverse();
		SizeType nl = left.size();

		SizeType offset = v.offset(0);
		SizeType npatches = ijPatch(GenIjPatchType::LEFT).size();

		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {

			SizeType igroup = ijPatch(GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = ijPatch(GenIjPatchType::RIGHT)[ipatch];

			SizeType sizeLeft = gengroupLeft_(igroup+1) - gengroupLeft_(igroup);
			SizeType sizeRight = gengroupRight_(jgroup+1) - gengroupRight_(jgroup);

			SizeType left_offset = gengroupLeft_(igroup);
			SizeType right_offset = gengroupRight_(jgroup);

			MatrixType& matrix = getMatrix(ipatch, sizeLeft, sizeRight);

			for (SizeType ileft=0; ileft < sizeLeft; ++ileft) {
				for (SizeType iright=0; iright < sizeRight; ++iright) {

					SizeType i = ileft + left_offset;
					SizeType j = iright + right_offset;

					SizeType ij = i + j * nl;

					assert(i < nl);
					assert(j < right.hamiltonian().row());

					assert(ij < permInverse.size());

					SizeType r = permInverse[ij];
					if (r < offset || r >= offset + v.effectiveSize(0))
						continue;

					matrix(ileft,iright) +=  v.slowAccess(r);
				}
			}
		}
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffsets<ComplexOrRealType>& v,
	                    SizeType targets)
	{
		err("useSvd doesn't yet work with VectorWithOffsets (sorry)\n");
	}

	MatrixType& getMatrix(SizeType ipatch, SizeType sizeLeft, SizeType sizeRight)
	{
		if (ipatch < allTargets_.size())
			return *(allTargets_[ipatch]);

		MatrixType* m = new MatrixType(sizeLeft, sizeRight);
		allTargets_.push_back(m);
		return *m;
	}

	void saveThisPatch(MatrixType& mAll,
	                   VectorRealType& eigs,
	                   const MatrixType& m,
	                   const MatrixType& vt,
	                   const VectorRealType& eigsOnePatch,
	                   SizeType ipatch)
	{
		GenIjPatchType& ijPatch = *(vectorOfijPatches_[0]);
		const MatrixType& mLeftOrRight = expandSys() ? m : vt;
		const GenGroupType& gengroup = expandSys() ? gengroupLeft_ : gengroupRight_;
		typename GenIjPatchType::LeftOrRightEnumType lOrR = (expandSys()) ?
		            GenIjPatchType::LEFT : GenIjPatchType::RIGHT;
		SizeType igroup = ijPatch(lOrR)[ipatch];
		SizeType offset = gengroup(igroup);
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
		return (params_.direction == ProgramGlobals::EXPAND_SYSTEM);
	}

	ProgressIndicatorType progress_;
	const ParamsType& params_;
	MatrixVectorType allTargets_;
	SparseMatrixType data_;
	bool debug_;
	bool verbose_;
	const LeftRightSuperType& lrs_;
	GenGroupType gengroupLeft_;
	GenGroupType gengroupRight_;
	VectorGenIjPatchType vectorOfijPatches_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

