/*
Copyright (c) 2012, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file InitKron.h
 *
 *
 */

#ifndef INIT_KRON_HEADER_H
#define INIT_KRON_HEADER_H

#include "ArrayOfMatStruct.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelType,typename ModelHelperType_>
class InitKron {

	static const bool KRON_USE_SYMMETRY = false;

public:

	typedef ModelHelperType_ ModelHelperType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef ArrayOfMatStruct<LeftRightSuperType> ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef typename ModelType::LinkProductStructType LinkProductStructType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

	InitKron(const ModelType& model,
	         const ModelHelperType& modelHelper)
	    : model_(model),
	      modelHelper_(modelHelper),
	      ijpatches_(modelHelper_.leftRightSuper(),modelHelper_.quantumNumber())
	{
		cacheSigns(modelHelper_.leftRightSuper().left().electronsVector());
		convertXcYcArrays();
		addHlAndHr();
	}

	~InitKron()
	{
		for (SizeType ic=0;ic<xc_.size();ic++) delete xc_[ic];
		for (SizeType ic=0;ic<yc_.size();ic++) delete yc_[ic];

	}

	bool useSymmetry() const { return KRON_USE_SYMMETRY; }

	bool loadBalance() const
	{
		return (model_.params().options.find("KronNoLoadBalance")
		        == PsimagLite::String::npos);
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

	const VectorSizeType& patch(typename GenIjPatchType::LeftOrRightEnumType i) const
	{
		return ijpatches_(i);
	}

	const LeftRightSuperType& lrs() const
	{
		return modelHelper_.leftRightSuper();
	}

	SizeType offset() const
	{
		SizeType m = modelHelper_.m();
		return modelHelper_.leftRightSuper().super().partition(m);
	}

	SizeType size() const { return modelHelper_.size(); }

	SizeType connections() const { return xc_.size(); }

	const ComplexOrRealType& value(SizeType i) const
	{
		assert(values_.size()>i);
		return values_[i];
	}

private:

	void addHlAndHr()
	{
		const RealType value = 1.0;
		const SparseMatrixType& aL = modelHelper_.leftRightSuper().left().hamiltonian();
		const SparseMatrixType& aR = modelHelper_.leftRightSuper().right().hamiltonian();
		identityL_.makeDiagonal(aL.row(), value);
		identityR_.makeDiagonal(aR.row(), value);
		std::pair<SizeType, SizeType> ops(0,0);
		std::pair<char, char> mods('n', 'n');
		LinkType link(0,0,0,value,0, ProgramGlobals::BOSON, ops, mods,1,value,0);
		addOneConnection(aL,identityR_,link);
		addOneConnection(identityL_,aR,link);
	}

	void convertXcYcArrays()
	{
		SizeType total = model_.getLinkProductStruct(modelHelper_);

		for (SizeType ix=0;ix<total;ix++) {
			SparseMatrixType const* A = 0;
			SparseMatrixType const* B = 0;

			LinkType link2 = model_.getConnection(&A,&B,ix,modelHelper_);
			if (link2.type==ProgramGlobals::ENVIRON_SYSTEM)  {
				LinkType link3 = link2;
				link3.type = ProgramGlobals::SYSTEM_ENVIRON;
				if (link3.fermionOrBoson == ProgramGlobals::FERMION)
					link3.value *= -1.0;
				addOneConnection(*B,*A,link3);
				continue;
			}

			addOneConnection(*A,*B,link2);
		}
	}

	void addOneConnection(const SparseMatrixType& A,
	                      const SparseMatrixType& B,
	                      const LinkType& link2)
	{
		SparseMatrixType Ahat;
		calculateAhat(Ahat, A, link2.value, link2.fermionOrBoson);
		values_.push_back(link2.value);
		ArrayOfMatStructType* x1 = new ArrayOfMatStructType(Ahat,
		                                                    ijpatches_,
		                                                    GenIjPatchType::LEFT);

		xc_.push_back(x1);

		ArrayOfMatStructType* y1 = new ArrayOfMatStructType(B,
		                                                    ijpatches_,
		                                                    GenIjPatchType::RIGHT);
		yc_.push_back(y1);
	}

	// Ahat(ia,ja) = (-1)^e_L(ia) A(ia,ja)*value
	void calculateAhat(SparseMatrixType& Ahat,
	                   const SparseMatrixType& A,
	                   ComplexOrRealType val,
	                   ProgramGlobals::FermionOrBosonEnum bosonOrFermion) const
	{
		Ahat = A;
		SizeType rows = Ahat.row();
		assert(signs_.size() == rows);
		SizeType counter = 0;
		for (SizeType i = 0; i < rows; ++i) {
			RealType sign = (bosonOrFermion == ProgramGlobals::FERMION &&
			                 signs_[i]) ? -1.0 : 1.0;
			for (int k = Ahat.getRowPtr(i); k < Ahat.getRowPtr(i+1); ++k) {
				ComplexOrRealType tmp = Ahat.getValue(k)*sign*val;
				Ahat.setValues(counter++, tmp);
			}
		}
	}

	void cacheSigns(const VectorSizeType& electrons)
	{
		signs_.resize(electrons.size(), false);
		for (SizeType i = 0; i < electrons.size(); ++i)
			signs_[i] = (electrons[i] & 1) ? true : false;
	}

	InitKron(const InitKron& other);

	InitKron& operator=(const InitKron& other);

	const ModelType& model_;
	const ModelHelperType& modelHelper_;
	GenIjPatchType  ijpatches_;
	SparseMatrixType identityL_;
	SparseMatrixType identityR_;
	VectorBoolType signs_;
	typename PsimagLite::Vector<ArrayOfMatStructType*>::Type xc_;
	typename PsimagLite::Vector<ArrayOfMatStructType*>::Type yc_;
	typename PsimagLite::Vector<ComplexOrRealType>::Type values_;

}; //class InitKron
} // namespace PsimagLite

/*@}*/

#endif // INIT_KRON_HEADER_H
