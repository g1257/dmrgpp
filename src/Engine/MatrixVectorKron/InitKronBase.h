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
/** \ingroup DMRG */
/*@{*/

/*! \file InitKronBase.h
 *
 *
 */
#ifndef INITKRON_BASE_H
#define INITKRON_BASE_H
#include "ProgramGlobals.h"
#include "ArrayOfMatStruct.h"
#include "Vector.h"
#include "Link.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class InitKronBase {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef Link<ComplexOrRealType> LinkType;
	typedef ArrayOfMatStruct<LeftRightSuperType> ArrayOfMatStructType;
	typedef typename LeftRightSuperType::BasisType BasisType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;

	enum WhatBasisEnum {OLD,  NEW};

	InitKronBase(const LeftRightSuperType& lrs,
	             SizeType m,
	             SizeType qn,
	             RealType denseSparseThreshold)
	    : mOld_(m),
	      mNew_(m),
	      denseSparseThreshold_(denseSparseThreshold),
	      ijpatchesOld_(lrs, qn),
	      ijpatchesNew_(&ijpatchesOld_),
	      wftMode_(false)
	{
		cacheSigns(signsNew_, lrs.left().electronsVector(BasisType::AFTER_TRANSFORM));
	}

	InitKronBase(const LeftRightSuperType& lrsOld,
	             SizeType mOld,
	             const LeftRightSuperType& lrsNew,
	             SizeType mNew,
	             SizeType qn,
	             RealType denseSparseThreshold)
	    : mOld_(mOld),
	      mNew_(mNew),
	      denseSparseThreshold_(denseSparseThreshold),
	      ijpatchesOld_(lrsOld, qn),
	      ijpatchesNew_(new GenIjPatchType(lrsNew, qn)),
	      wftMode_(true)
	{
		cacheSigns(signsNew_, lrsNew.left().electronsVector(BasisType::AFTER_TRANSFORM));
	}

	~InitKronBase()
	{
		for (SizeType ic=0;ic<xc_.size();ic++) delete xc_[ic];
		for (SizeType ic=0;ic<yc_.size();ic++) delete yc_[ic];
		if (wftMode_) {
			delete ijpatchesNew_;
			ijpatchesNew_ = 0;
		}
	}

	const LeftRightSuperType& lrs(WhatBasisEnum what) const
	{
		return (what == OLD) ? ijpatchesOld_.lrs() : ijpatchesNew_->lrs();
	}

	const VectorSizeType& patch(WhatBasisEnum what,
	                            typename GenIjPatchType::LeftOrRightEnumType i) const
	{
		return (what == OLD) ? ijpatchesOld_(i) :
		                       ijpatchesNew_->operator()(i);
	}

	SizeType offset(WhatBasisEnum what) const
	{
		return (what == OLD) ? ijpatchesOld_.lrs().super().partition(mOld_) :
		                       ijpatchesNew_->lrs().super().partition(mNew_);
	}

	const ArrayOfMatStructType& xc(SizeType ic) const
	{
		assert(ic < xc_.size());
		assert(xc_[ic]);
		return *xc_[ic];
	}

	const ArrayOfMatStructType& yc(SizeType ic) const
	{
		assert(ic < yc_.size());
		assert(yc_[ic]);
		return *yc_[ic];
	}

	const ComplexOrRealType& value(SizeType i) const
	{
		assert(values_.size()>i);
		return values_[i];
	}

	SizeType connections() const { return xc_.size(); }

	SizeType size(WhatBasisEnum what) const
	{
		return (what == OLD) ? sizeInternal(ijpatchesOld_, mOld_) :
		                       sizeInternal(*ijpatchesNew_, mNew_);
	}


	const VectorSizeType& weightsOfPatchesNew() const
	{
		return weightsOfPatches_;
	}


	void computeOffsets(VectorSizeType& offsetForPatches,
	                    WhatBasisEnum what)
	{
		assert(patch(what, GenIjPatchType::LEFT).size() ==
		       patch(what, GenIjPatchType::RIGHT).size());

		SizeType npatch = patch(what, GenIjPatchType::LEFT).size();
		SizeType sum = 0;
		const BasisType& left = lrs(what).left();
		const BasisType& right = lrs(what).right();

		for( SizeType ipatch=0; ipatch < npatch; ipatch++) {
			offsetForPatches[ipatch] = sum;

			SizeType igroup = patch(what, GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = patch(what, GenIjPatchType::RIGHT)[ipatch];

			assert(left.partition(igroup+1) >= left.partition(igroup));
			SizeType sizeLeft =  left.partition(igroup+1) - left.partition(igroup);
			assert(right.partition(jgroup+1) >= right.partition(jgroup));
			SizeType sizeRight = right.partition(jgroup+1) - right.partition(jgroup);

			sum +=  sizeLeft*sizeRight;
		}
	}

	SizeType numberOfPatches(WhatBasisEnum what) const
	{
		assert(patch(what, GenIjPatchType::LEFT).size() ==
		       patch(what, GenIjPatchType::RIGHT).size());
		return patch(what, GenIjPatchType::LEFT).size();
	}
protected:

	void addOneConnection(const SparseMatrixType& A,
	                      const SparseMatrixType& B,
	                      const LinkType& link2)
	{
		SparseMatrixType Ahat;
		calculateAhat(Ahat, A, link2.value, link2.fermionOrBoson);
		values_.push_back(link2.value);
		ArrayOfMatStructType* x1 = new ArrayOfMatStructType(Ahat,
		                                                    ijpatchesOld_,
		                                                    *ijpatchesNew_,
		                                                    GenIjPatchType::LEFT,
		                                                    denseSparseThreshold_);

		xc_.push_back(x1);

		ArrayOfMatStructType* y1 = new ArrayOfMatStructType(B,
		                                                    ijpatchesOld_,
		                                                    *ijpatchesNew_,
		                                                    GenIjPatchType::RIGHT,
		                                                    denseSparseThreshold_);
		yc_.push_back(y1);
	}

	// -------------------------------------------
	// setup vstart(:) for beginning of each patch
	// -------------------------------------------
	void setUpVstart(VectorSizeType& vstart, WhatBasisEnum what)
	{
		SizeType npatches = patch(what, GenIjPatchType::LEFT).size();
		assert(npatches > 0);
		SizeType ip = 0;
		VectorSizeType weights(npatches, 0);
		const BasisType& left = lrs(what).left();
		const BasisType& right = lrs(what).right();

		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {
			vstart[ipatch] = ip;

			SizeType igroup = patch(what, GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = patch(what, GenIjPatchType::RIGHT)[ipatch];

			assert(left.partition(igroup+1) >= left.partition(igroup));
			SizeType sizeLeft =  left.partition(igroup+1) - left.partition(igroup);

			assert(right.partition(jgroup+1) >= right.partition(jgroup));
			SizeType sizeRight = right.partition(jgroup+1) - right.partition(jgroup);

			assert(1 <= sizeLeft);
			assert(1 <= sizeRight);

			weights[ipatch] = sizeLeft * sizeRight * (sizeLeft + sizeRight);

			ip += sizeLeft * sizeRight;
		}

		vstart[npatches] = ip;

		if (what == NEW)
			setAndFixWeights(weights);
	}

	// -------------------
	// copy xout(:) to vout(:)
	// -------------------
	void copyOut(VectorType& vout,
	             const VectorType& xout,
	             const VectorSizeType& vstart) const
	{
		const VectorSizeType& permInverse = lrs(NEW).super().permutationInverse();
		SizeType offset1 = offset(NEW);
		SizeType nl = lrs(NEW).left().hamiltonian().rows();
		SizeType npatches = patch(NEW, GenIjPatchType::LEFT).size();
		const BasisType& left = lrs(NEW).left();
		const BasisType& right = lrs(NEW).right();

		for( SizeType ipatch=0; ipatch < npatches; ++ipatch) {

			SizeType igroup = patch(NEW, GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = patch(NEW, GenIjPatchType::RIGHT)[ipatch];

			assert(left.partition(igroup+1) >= left.partition(igroup));
			SizeType sizeLeft =  left.partition(igroup+1) - left.partition(igroup);

			assert(right.partition(jgroup+1) >= right.partition(jgroup));
			SizeType sizeRight = right.partition(jgroup+1) - right.partition(jgroup);

			SizeType left_offset = left.partition(igroup);
			SizeType right_offset = right.partition(jgroup);

			for (SizeType ileft=0; ileft < sizeLeft; ++ileft) {
				for (SizeType iright=0; iright < sizeRight; ++iright) {

					SizeType i = ileft + left_offset;
					SizeType j = iright + right_offset;

					assert(i < nl);
					assert(j < lrs(NEW).right().hamiltonian().rows());

					assert(i + j*nl < permInverse.size());

					SizeType r = permInverse[i + j*nl];
					assert( !(  (r < offset1) || (r >= (offset1 + size(NEW))) ) );

					SizeType ip = vstart[ipatch] + (iright + ileft * sizeRight);
					assert(ip < xout.size());

					assert(r >= offset1 && ((r - offset1) < vout.size()) );
					vout[r-offset1] = xout[ip];
				}
			}
		}
	}

private:

	void setAndFixWeights(const VectorSizeType& weights)
	{
		long unsigned int max = *(std::max_element(weights.begin(), weights.end()));
		max >>= 31;
		SizeType bits = 1 + PsimagLite::log2Integer(max);
		SizeType npatches = weights.size();
		weightsOfPatches_.resize(npatches);
		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {
			long unsigned int tmp = (weights[ipatch] >> bits);
			weightsOfPatches_[ipatch] = (max == 0) ? weights[ipatch] : tmp;
		}
	}

	static SizeType sizeInternal(const GenIjPatchType& ijpatches,
	                             SizeType m)
	{
		assert(ijpatches.lrs().super().partition(m + 1) >=
		       ijpatches.lrs().super().partition(m));
		return ijpatches.lrs().super().partition(m + 1) -
		        ijpatches.lrs().super().partition(m);
	}

	static void cacheSigns(VectorBoolType& signs, const VectorSizeType& electrons)
	{
		signs.resize(electrons.size(), false);
		for (SizeType i = 0; i < electrons.size(); ++i)
			signs[i] = (electrons[i] & 1) ? true : false;
	}

	// Ahat(ia,ja) = (-1)^e_L(ia) A(ia,ja)*value
	void calculateAhat(SparseMatrixType& Ahat,
	                   const SparseMatrixType& A,
	                   ComplexOrRealType val,
	                   ProgramGlobals::FermionOrBosonEnum bosonOrFermion) const
	{
		Ahat = A;
		SizeType rows = Ahat.rows();
		assert(signsNew_.size() >= rows);
		SizeType counter = 0;
		for (SizeType i = 0; i < rows; ++i) {
			RealType sign = (bosonOrFermion == ProgramGlobals::FERMION &&
			                 signsNew_[i]) ? -1.0 : 1.0;
			for (int k = Ahat.getRowPtr(i); k < Ahat.getRowPtr(i+1); ++k) {
				ComplexOrRealType tmp = Ahat.getValue(k)*sign*val;
				Ahat.setValues(counter++, tmp);
			}
		}
	}

	InitKronBase(const InitKronBase&);

	InitKronBase& operator=(const InitKronBase&);

	SizeType mOld_;
	SizeType mNew_;
	const RealType denseSparseThreshold_;
	GenIjPatchType ijpatchesOld_;
	GenIjPatchType* ijpatchesNew_;
	VectorSizeType weightsOfPatches_;
	VectorArrayOfMatStructType xc_;
	VectorArrayOfMatStructType yc_;
	VectorType values_;
	VectorBoolType signsNew_;
	bool wftMode_;
};
} // namespace Dmrg

#endif // INITKRON_BASE_H
