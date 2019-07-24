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

/*! \file InitKronHamiltonian.h
 *
 *
 */
#ifndef INITKRON_HAMILTONIAN_H
#define INITKRON_HAMILTONIAN_H
#include "ProgramGlobals.h"
#include "InitKronBase.h"
#include "Vector.h"
#include "Profiling.h"

namespace Dmrg {

template<typename ModelType_>
class InitKronHamiltonian : public InitKronBase<typename ModelType_::LeftRightSuperType> {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef ModelType_ ModelType;
	typedef typename ModelType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename ModelHelperType::OperatorStorageType OperatorStorageType;
	typedef typename LeftRightSuperType::BasisType BasisType;
	typedef InitKronBase<LeftRightSuperType> BaseType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BaseType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;

	InitKronHamiltonian(const ModelType& model,
	                    const HamiltonianConnectionType& hc)
	    : BaseType(hc.modelHelper().leftRightSuper(),
	               hc.modelHelper().m(),
	               hc.modelHelper().quantumNumber(),
	               model.params().denseSparseThreshold,
	               model.params().options.find("KronNoUseLowerPart") == PsimagLite::String::npos
	               && model.params().options.find("BatchedGemm") == PsimagLite::String::npos),
	      model_(model),
	      hc_(hc),
	      vstart_(BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size() + 1),
	      offsetForPatches_(BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size() + 1)
	{
		addHlAndHr();

		{
			PsimagLite::Profiling profiling("convertXcYcArrays", std::cout);

			convertXcYcArrays();
		}

		BaseType::setUpVstart(vstart_, BaseType::NEW);
		assert(vstart_.size() > 0);
		SizeType nsize = vstart_[vstart_.size() - 1];
		assert(nsize > 0);
		yin_.resize(nsize, 0.0);
		xout_.resize(nsize, 0.0);
		BaseType::computeOffsets(offsetForPatches_, BaseType::NEW);
	}

	bool isWft() const {return false; }

	bool loadBalance() const
	{
		return (model_.params().options.find("KronLoadBalance") != PsimagLite::String::npos);
	}

	// -------------------
	// copy vin(:) to yin(:)
	// -------------------
	void copyIn(const VectorType& vout,
	            const VectorType& vin)
	{
		VectorType& xout = xout_;
		VectorType& yin = yin_;

		const VectorSizeType& permInverse = BaseType::lrs(BaseType::NEW).super().permutationInverse();
		const SparseMatrixType& leftH = BaseType::lrs(BaseType::NEW).left().hamiltonian();
		SizeType nl = leftH.rows();

		SizeType offset = BaseType::offset(BaseType::NEW);
		SizeType npatches = BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT).size();
		const BasisType& left = BaseType::lrs(BaseType::NEW).left();
		const BasisType& right = BaseType::lrs(BaseType::NEW).right();

		for (SizeType ipatch=0; ipatch < npatches; ++ipatch) {

			SizeType igroup = BaseType::patch(BaseType::NEW, GenIjPatchType::LEFT)[ipatch];
			SizeType jgroup = BaseType::patch(BaseType::NEW, GenIjPatchType::RIGHT)[ipatch];

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
					assert(j < BaseType::lrs(BaseType::NEW).right().hamiltonian().rows());

					assert(ij < permInverse.size());

					SizeType r = permInverse[ ij ];
					assert(!((r < offset) || (r >= (offset + BaseType::size(BaseType::NEW)))));

					SizeType ip = vstart_[ipatch] + (iright + ileft * sizeRight);
					assert(ip < yin.size());

					assert( (r >= offset) && ((r-offset) < vin.size()) );
					yin[ip] = vin[r-offset];
					xout[ip] = vout[r-offset];
				}
			}
		}
	}

	// -------------------
	// copy xout(:) to vout(:)
	// -------------------
	void copyOut(VectorType& vout) const
	{
		BaseType::copyOut(vout, xout_, vstart_);
	}

	const VectorType& yin() const { return yin_; }

	VectorType& xout() { return xout_; }

	const SizeType& offsetForPatches(typename BaseType::WhatBasisEnum,
	                                 SizeType ind) const
	{
		assert(ind < offsetForPatches_.size());
		return  offsetForPatches_[ind];
	}

	bool batchedGemm() const
	{
		return (model_.params().options.find("BatchedGemm") != PsimagLite::String::npos);
	}

private:

	void addHlAndHr()
	{
		const RealType value = 1.0;
		const OperatorStorageType& aL = hc_.modelHelper().leftRightSuper().left().hamiltonian();
		const OperatorStorageType& aR = hc_.modelHelper().leftRightSuper().right().hamiltonian();
		identityL_.makeDiagonal(aL.rows(), value);
		identityR_.makeDiagonal(aR.rows(), value);
		std::pair<SizeType, SizeType> ops(0,0);
		std::pair<char, char> mods('n', 'n');
		LinkType link(0,
		              0,
		              ProgramGlobals::ConnectionEnum::SYSTEM_SYSTEM,
		              value,
		              ProgramGlobals::FermionOrBosonEnum::BOSON,
		              ops,
		              mods,
		              1,
		              value,
		              0);
		BaseType::addOneConnection(aL,identityR_,link);
		BaseType::addOneConnection(identityL_,aR,link);
	}

	void convertXcYcArrays()
	{
		SizeType total = hc_.tasks();

		for (SizeType ix=0;ix<total;ix++) {
			OperatorStorageType const* A = 0;
			OperatorStorageType const* B = 0;

			LinkType link2 = hc_.getKron(&A, &B, ix);
			if (link2.type==ProgramGlobals::ConnectionEnum::ENVIRON_SYSTEM)  {
				LinkType link3 = link2;
				link3.type = ProgramGlobals::ConnectionEnum::SYSTEM_ENVIRON;
				if (link3.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION)
					link3.value *= -1.0;

				assert(A);
				assert(B);
				BaseType::addOneConnection(*B,*A,link3);
				continue;
			}

			assert(A);
			assert(B);
			BaseType::addOneConnection(*A,*B,link2);
		}
	}

	InitKronHamiltonian(const InitKronHamiltonian&);

	InitKronHamiltonian& operator=(const InitKronHamiltonian&);

	const ModelType& model_;
	const HamiltonianConnectionType& hc_;
	OperatorStorageType identityL_;
	OperatorStorageType identityR_;
	VectorSizeType vstart_;
	VectorType yin_;
	VectorType xout_;
	VectorSizeType offsetForPatches_;
};
} // namespace Dmrg

#endif // INITKRON_HAMILTONIAN_H
