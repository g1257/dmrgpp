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
	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef PsimagLite::ProgressIndicator ProgressIndicatorType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename BaseType::Params ParamsType;
	typedef GenIjPatch<LeftRightSuperType> GenIjPatchType;
	typedef typename GenIjPatchType::VectorSizeType VectorSizeType;
	typedef typename GenIjPatchType::GenGroupType GenGroupType;

	enum {EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM };

public:

	DensityMatrixSvd(const TargettingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    :
	      progress_("DensityMatrixSvd"),
	      debug_(p.debug),
	      verbose_(p.verbose)
	{
		{
			PsimagLite::OstringStream msg;
			msg<<"Init partition for all targets";
			progress_.printline(msg,std::cout);
		}

		const BasisWithOperatorsType& left = lrs.left();
		const BasisWithOperatorsType& right = lrs.right();
		SizeType oneOrZero = (target.includeGroundStage()) ? 1 : 0;
		SizeType targets = oneOrZero + target.size(); // Number of targets;
		SizeType freeSize = (p.direction == ProgramGlobals::EXPAND_SYSTEM) ?
		            left.size() : right.size();
		SizeType summedSize = (p.direction == ProgramGlobals::EXPAND_SYSTEM) ?
		            right.size() : left.size();
		allTargets_.resize(targets*summedSize, freeSize);
		allTargets_.setTo(0.0);
		GenGroupType gengroupLeft(left);
		GenGroupType gengroupRight(right);
		for (SizeType x = 0; x < targets; ++x)
			addThisTarget(x, target, gengroupLeft, gengroupRight, lrs, targets);

		{
			PsimagLite::OstringStream msg;
			msg<<"Done with init partition";
			progress_.printline(msg,std::cout);
		}
	}

	virtual SparseMatrixType& operator()()
	{
		return data_;
	}

	virtual SizeType rank() { return allTargets_.rows(); }

	void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz)
	{
		SizeType freeSize = allTargets_.cols();
		MatrixType vt(freeSize, freeSize);
		eigs.resize(freeSize);
		svd('A', allTargets_, eigs, vt);
		fullMatrixToCrsMatrix(data_, allTargets_);
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
	                   const GenGroupType& genGroupLeft,
	                   const GenGroupType& genGroupRight,
	                   const LeftRightSuperType& lrs,
	                   SizeType targets)

	{
		SizeType x2 = (target.includeGroundStage() && x > 0 ) ? x - 1 : x;

		const VectorWithOffsetType& v = (target.includeGroundStage() && x == 0) ?
		            target.gs() : target(x2);

		addThisTarget2(x, v, genGroupLeft, genGroupRight, lrs, targets);
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffset<ComplexOrRealType>& v,
	                    const GenGroupType& genGroupLeft,
	                    const GenGroupType& genGroupRight,
	                    const LeftRightSuperType& lrs,
	                    SizeType targets)
	{
		const BasisType& super = lrs.super();
		const BasisWithOperatorsType& left = lrs.left();
		const BasisWithOperatorsType& right = lrs.right();

		SizeType m = v.sector(0);
		int state = super.partition(m);
		SizeType qn = super.qn(state);
		GenIjPatchType ijPatch(lrs, qn);

		const VectorSizeType& permInverse = super.permutationInverse();
		SizeType nl = left.size();

		SizeType offset = v.offset(0);
		SizeType npatches = ijPatch(GenIjPatchType::LEFT).size();

		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {

			SizeType igroup = ijPatch(GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = ijPatch(GenIjPatchType::RIGHT)[ipatch];

			SizeType sizeLeft = genGroupLeft(igroup+1) - genGroupLeft(igroup);
			SizeType sizeRight = genGroupRight(jgroup+1) - genGroupRight(jgroup);

			SizeType left_offset = genGroupLeft(igroup);
			SizeType right_offset = genGroupRight(jgroup);

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

					allTargets_(i,j) =  v.slowAccess(r);

				}
			}
		}
	}

	void addThisTarget2(SizeType x,
	                    const VectorWithOffsets<ComplexOrRealType>& v,
	                    const BasisWithOperatorsType& left,
	                    const BasisWithOperatorsType& right,
	                    const BasisType& super,
	                    SizeType targets)
	{
		err("useSvd doesn't yet work with VectorWithOffsets (sorry)\n");
	}

	ProgressIndicatorType progress_;
	MatrixType allTargets_;
	SparseMatrixType data_;
	bool debug_;
	bool verbose_;
}; // class DensityMatrixSvd

} // namespace Dmrg

#endif

