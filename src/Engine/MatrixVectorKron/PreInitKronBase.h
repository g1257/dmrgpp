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

/*! \file PreInitKronBase.h
 *
 *
 */
#ifndef PREINITKRON_BASE_H
#define PREINITKRON_BASE_H
#include "ProgramGlobals.h"
#include "ArrayOfMatStruct.h"
#include "Vector.h"
#include "Link.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class PreInitKronBase {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;

public:

	typedef typename LeftRightSuperType::SparseMatrixType SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef Link<ComplexOrRealType> LinkType;
	typedef ArrayOfMatStruct<LeftRightSuperType> ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::GenIjPatchType GenIjPatchType;
	typedef typename PsimagLite::Vector<ArrayOfMatStructType*>::Type VectorArrayOfMatStructType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename ArrayOfMatStructType::VectorSizeType VectorSizeType;

	PreInitKronBase(const LeftRightSuperType& lrs,
	                SizeType m,
	                SizeType qn,
	                RealType denseSparseThreshold)
	    : m_(m),
	      denseSparseThreshold_(denseSparseThreshold),
	      ijpatches_(lrs, qn)
	{
		cacheSigns(lrs.left().electronsVector());
	}

	~PreInitKronBase()
	{
		for (SizeType ic=0;ic<xc_.size();ic++) delete xc_[ic];
		for (SizeType ic=0;ic<yc_.size();ic++) delete yc_[ic];
	}

	const LeftRightSuperType& lrs() const
	{
		return ijpatches_.lrs();
	}

	const VectorSizeType& patch(typename GenIjPatchType::LeftOrRightEnumType i) const
	{
		return ijpatches_(i);
	}

	SizeType offset() const
	{
		return ijpatches_.lrs().super().partition(m_);
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

	SizeType size() const
	{
		assert(ijpatches_.lrs().super().partition(m_ + 1) >=
		       ijpatches_.lrs().super().partition(m_));
		return ijpatches_.lrs().super().partition(m_ + 1) -
		       ijpatches_.lrs().super().partition(m_);
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
		                                                    ijpatches_,
		                                                    GenIjPatchType::LEFT,
		                                                    denseSparseThreshold_);

		xc_.push_back(x1);

		ArrayOfMatStructType* y1 = new ArrayOfMatStructType(B,
		                                                    ijpatches_,
		                                                    GenIjPatchType::RIGHT,
		                                                    denseSparseThreshold_);
		yc_.push_back(y1);
	}

private:

	// Ahat(ia,ja) = (-1)^e_L(ia) A(ia,ja)*value
	void calculateAhat(SparseMatrixType& Ahat,
	                   const SparseMatrixType& A,
	                   ComplexOrRealType val,
	                   ProgramGlobals::FermionOrBosonEnum bosonOrFermion) const
	{
		Ahat = A;
		SizeType rows = Ahat.rows();
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

	PreInitKronBase(const PreInitKronBase&);

	PreInitKronBase& operator=(const PreInitKronBase&);

	SizeType m_;
	const RealType& denseSparseThreshold_;
	GenIjPatchType ijpatches_;
	VectorArrayOfMatStructType xc_;
	VectorArrayOfMatStructType yc_;
	VectorType values_;
	VectorBoolType signs_;
};
} // namespace Dmrg

#endif // PREINITKRON_BASE_H
