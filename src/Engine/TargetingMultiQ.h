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
	typedef typename BaseType::VectorVectorVectorType VectorVectorVectorType;

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
	      weight_(vqn_.size())
	{}

	virtual bool includeGroundStage() const { return false; }

	virtual void set(VectorVectorVectorType& inV,
	                 const VectorSizeType& sectors,
	                 const BasisType& basis)
	{
		const SizeType n = sectors.size();
		if (inV.size() != 1)
			err("Expected exactly one vector instead of " + ttos(inV.size()) +
			    + ", that is,  just the ground state\n");

		if (inV[0].size() != n)
			err("Expected inV[0].size == " + ttos(sectors) + " but found "
			    + ttos(inV[0].size()) + " instead\n");

		if (this->common().aoe().targetVectors().size() != n)
			this->common().aoe().targetVectorsResize(n);

		for (SizeType i = 0; i < n; ++i) {
			SizeType j = sectors[i];
			VectorSizeType weights(basis.partition() - 1);

			weights[j] = basis.partition(j + 1) - basis.partition(j);
			VectorWithOffsetType vwo(weights, basis);
			vwo.setDataInSector(inV[0][i], j);
			VectorWithOffsetType& handle =
			        const_cast<VectorWithOffsetType&>(this->common().aoe().targetVectors()[i]);
			handle = vwo;
			weight_[i] = 1.0/n;
		}
	}

	void initialGuess(typename PsimagLite::Vector<VectorType>::Type& initialVector,
	                  const VectorSizeType& block,
	                  bool noguess,
	                  VectorSizeType& weights,
	                  const BasisType& basis) const
	{
		VectorSizeType sectors;
		findSectors(sectors, weights);
		const SizeType n = sectors.size();
		initialVector.resize(n);
		for (SizeType i = 0; i < n; ++i) {
			VectorSizeType weights2(weights.size());
			weights2[sectors[i]] = weights[sectors[i]];
			VectorWithOffsetType vwo(weights2, basis);
			this->common().initialGuess(vwo,
			                            this->common().aoe().targetVectors()[i],
			                            block,
			                            noguess);
			vwo.extract(initialVector[i], sectors[i]);
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

	void evolve(const VectorRealType&,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType&,
	            SizeType)
	{
		if (this->common().aoe().targetVectors().size() != weight_.size())
			return;

		const bool doBorderIfBorder = true;
		this->common().cocoon(block1, direction, doBorderIfBorder);
	}

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

	void findSectors(VectorSizeType& sectors, const VectorSizeType& weights) const
	{
		const SizeType m = weights.size();
		for (SizeType i = 0; i < m; ++i)
			if (weights[i] > 0) sectors.push_back(i);
	}

	const VectorQnType& vqn_;
	TargetParamsType tstStruct_;
	PsimagLite::ProgressIndicator progress_;
	VectorRealType weight_;
};     //class TargetingMultiQ
} // namespace Dmrg
/*@}*/
#endif

