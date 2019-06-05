/*
Copyright (c) 2009-2019, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

/*! \file TargetingMultiQ.h
 *
 *  targets the ground state
 *
 */

#ifndef TARGETING_MULTI_Q
#define TARGETING_MULTI_Q
#include <iostream>
#include "TargetParamsGroundState.h"
#include "ApplyOperatorLocal.h"
#include <stdexcept>
#include "PsimagLite.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "TargetQuantumElectrons.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingMultiQ : public TargetingBase<LanczosSolverType_, VectorWithOffsetType_> {

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType, VectorWithOffsetType_> BaseType;
	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef VectorType TargetVectorType;
	typedef TargetParamsGroundState<ModelType> TargetParamsType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename BasisType::QnType QnType;
	typedef typename QnType::VectorQnType VectorQnType;
	typedef TargetQuantumElectrons<RealType, QnType> TargetQuantumElectronsType;

	TargetingMultiQ(const LeftRightSuperType& lrs,
	                const ModelType& model,
	                const WaveFunctionTransfType& wft,
	                const VectorQnType& vqn,
	                InputValidatorType&,
	                PsimagLite::String targeting)
	    : BaseType(lrs, model, wft, 0),
	      vqn_(vqn),
	      tstStruct_(targeting),
	      progress_(targeting),
	      gsWeight_(tstStruct_.gsWeight())
	{}

	virtual bool includeGroundStage() const { return false; }

	virtual void set(typename PsimagLite::Vector<VectorType>::Type& v,
	                 const VectorSizeType& sectors,
	                 const BasisType& basis)
	{
		const SizeType n = sectors.size();
		if (this->common().aoe().targetVectors().size() == 0)
			this->common().aoe().targetVectorsResize(n);
		else
			if (this->common().aoe().targetVectors().size() != n)
				err("TargetingMultiQ: Wrong number of targets (FATAL ERROR)\n");

		VectorSizeType weights(basis.partition() - 1);
		for (SizeType i = 0; i < n; ++i) {
			SizeType j = sectors[i];
			weights[j] = basis.partition(j + 1) - basis.partition(j);
		}

		const RealType factor = 1 - gsWeight_;
		for (SizeType i = 0; i < n; ++i) {
			VectorWithOffsetType vwo(weights, basis);
			vwo.setDataInSector(v[i], i);
			VectorWithOffsetType& handle =
			        const_cast<VectorWithOffsetType&>(this->common().aoe().targetVectors()[i]);
			handle = vwo;
			weight_[i] = factor/n;
		}
	}

	SizeType sites() const { return 0; }

	SizeType targets() const { return weight_.size(); }

	RealType weight(SizeType i) const
	{
		assert(i < weight_.size());
		return  weight_[i];
	}

	RealType gsWeight() const
	{
		throw PsimagLite::RuntimeError("gsWeight should not be called by this target\n");
	}

	SizeType size() const { return weight_.size(); }

	void evolve(RealType,
	            ProgramGlobals::DirectionEnum,
	            const BlockType&,
	            const BlockType&,
	            SizeType)
	{}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
	}

	void read(typename TargetingCommonType::IoInputType& io,
	          PsimagLite::String prefix)
	{
		this->common().read(io, prefix);
	}

private:

	const VectorQnType& vqn_;
	TargetParamsType tstStruct_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	VectorRealType weight_;
};     //class TargetingMultiQ
} // namespace Dmrg
/*@}*/
#endif

