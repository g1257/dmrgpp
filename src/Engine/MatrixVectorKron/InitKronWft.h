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

/*! \file InitKronWft.h
 *
 *
 */
#ifndef INITKRON_WFT_H
#define INITKRON_WFT_H
#include "ProgramGlobals.h"
#include "InitKronBase.h"
#include "Vector.h"

namespace Dmrg {

template<typename LeftRightSuperType_, typename WftOptionsType, typename DmrgWaveStructType>
class InitKronWft : public InitKronBase<LeftRightSuperType_> {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef InitKronBase<LeftRightSuperType> BaseType;
	typedef typename LeftRightSuperType::BasisType BasisType;
	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename BaseType::LinkType LinkType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BaseType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;
	typedef typename LinkType::PairSizeType PairSizetype;
	typedef typename LinkType::PairCharType PairCharType;

	InitKronWft(const LeftRightSuperType& lrsNew,
	            SizeType mNew,
	            SizeType qn,
	            const WftOptionsType& wftOptions,
	            const DmrgWaveStructType& dmrgWaveStruct,
	            SizeType mOld)
	    : BaseType(dmrgWaveStruct.lrs,
	               mOld,
	               lrsNew,
	               mNew,
	               qn,
	               wftOptions.denseSparseThreshold),
	      wftOptions_(wftOptions),
	      vstartNew_(BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size() + 1),
	      vstartOld_(BaseType::patch(BaseType::OLD, GenIjPatchType::LEFT).size() + 1),
	      offsetForPatchesNew_(BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size()),
	      offsetForPatchesOld_(BaseType::patch(BaseType::OLD, GenIjPatchType::LEFT).size())
	{
		BaseType::setUpVstart(vstartNew_, BaseType::NEW);
		assert(vstartNew_.size() > 0);
		SizeType nsizeNew = vstartNew_[vstartNew_.size() - 1];
		assert(nsizeNew > 0);
		xout_.resize(nsizeNew, 0.0);

		BaseType::setUpVstart(vstartOld_, BaseType::OLD);
		assert(vstartOld_.size() > 0);
		SizeType nsizeOld = vstartOld_[vstartOld_.size() - 1];
		assert(nsizeOld > 0);
		yin_.resize(nsizeOld, 0.0);

		SparseMatrixType we;
		dmrgWaveStruct.we.toSparse(we);
		const PairCharType nn('N', 'N');
		const PairSizetype dummy(0, 0);
		const ComplexOrRealType value = 1.0;

		LinkType link(0,
		              0,
		              ProgramGlobals::SYSTEM_ENVIRON,
		              value,
		              1,
		              ProgramGlobals::BOSON,
		              dummy,
		              nn,
		              1,
		              1,
		              0);

		SparseMatrixType ws;
		dmrgWaveStruct.ws.toSparse(ws);
		SparseMatrixType weT;
		transposeConjugate(weT,we);
		SparseMatrixType wsT;
		transposeConjugate(wsT,ws);
		if (wftOptions_.dir == ProgramGlobals::EXPAND_SYSTEM) {

			BaseType::addOneConnection(wsT, we, link);
		} else {

			BaseType::addOneConnection(ws, we, link);
		}

		BaseType::computeOffsets(offsetForPatchesNew_, BaseType::NEW);
		BaseType::computeOffsets(offsetForPatchesOld_, BaseType::OLD);
	}

	bool isWft() const {return true; }

	bool loadBalance() const
	{
		return wftOptions_.kronLoadBalance;
	}

	// -------------------
	// copy vin(:) to yin(:)
	// -------------------
	void copyIn(const VectorType& vout,
	            const VectorType& vin)
	{
		copyIn(xout_, vout, vstartNew_, BaseType::NEW);
		copyIn(yin_, vin, vstartOld_, BaseType::OLD);
	}

	// -------------------
	// copy xout(:) to vout(:)
	// -------------------
	void copyOut(VectorType& vout) const
	{
		BaseType::copyOut(vout, xout_, vstartNew_);
	}

	const VectorType& yin() const { return yin_; }

	VectorType& xout() { return xout_; }

	SizeType offsetForPatches(typename BaseType::WhatBasisEnum what) const
	{
		return (what == BaseType::NEW) ? offsetForPatchesNew_.size() :
		                                 offsetForPatchesOld_.size();
	}

	const SizeType& offsetForPatches(typename BaseType::WhatBasisEnum what,
	                                       SizeType ind) const
	{
		assert(ind < offsetForPatches(what));
		return (what == BaseType::NEW) ? offsetForPatchesNew_[ind] :
		                                 offsetForPatchesOld_[ind];
	}

private:

	void copyIn(VectorType& x,
	            const VectorType& v,
	            const VectorSizeType& vstart,
	            typename BaseType::WhatBasisEnum what) const
	{
		const VectorSizeType& permInverse = BaseType::lrs(what).super().permutationInverse();
		const SparseMatrixType& leftH = BaseType::lrs(what).left().hamiltonian();
		SizeType nl = leftH.rows();

		SizeType offset1 = BaseType::offset(what);
		SizeType npatches = BaseType::patch(what, GenIjPatchType::LEFT).size();
		const BasisType& left = BaseType::lrs(what).left();
		const BasisType& right = BaseType::lrs(what).right();

		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {

			SizeType igroup = BaseType::patch(what, GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = BaseType::patch(what, GenIjPatchType::RIGHT)[ipatch];

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

					SizeType ij = i + j * nl;

					assert(i < nl);
					assert(j < BaseType::lrs(what).right().hamiltonian().rows());

					assert(ij < permInverse.size());

					SizeType r = permInverse[ ij ];
					assert(ipatch < vstart.size());
					SizeType ip = vstart[ipatch] + (iright + ileft * sizeRight);
					assert(ip < x.size());
					assert(r >= offset1);
					r -= offset1;
					assert(r < v.size());

					x[ip] = v[r];
				}
			}
		}
	}

	InitKronWft(const InitKronWft&);

	InitKronWft& operator=(const InitKronWft&);

	const WftOptionsType& wftOptions_;
	VectorSizeType vstartNew_;
	VectorSizeType vstartOld_;
	VectorType yin_;
	VectorType xout_;
	VectorSizeType offsetForPatchesNew_;
	VectorSizeType offsetForPatchesOld_;
};
} // namespace Dmrg

#endif // INITKRON_WFT_H
