/*
Copyright (c) 2009-2011, UT-Battelle, LLC
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

/** \ingroup DMRG */
/*@{*/

/*! \file TargetingCorrectionVector.h
 *
 * Implements the targeting required by
 * the correction targeting method
 *
 */

#ifndef TARGETING_CORRECTION_VECTOR_H
#define TARGETING_CORRECTION_VECTOR_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TargetParamsCorrectionVector.h"
#include "VectorWithOffsets.h"
#include "CorrectionVectorFunction.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "FreqEnum.h"
#include "NoPthreadsNg.h"
#include "CorrectionVectorSkeleton.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingCorrectionVector : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;

public:

	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef TargetParamsCorrectionVector<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef typename TargetingCommonType::TimeSerializerType TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename LanczosSolverType::MatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,
	TargetParamsType> CorrectionVectorFunctionType;
	typedef ParallelTriDiag<ModelType,
	LanczosSolverType,
	VectorWithOffsetType> ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef CorrectionVectorSkeleton<LanczosSolverType,
	VectorWithOffsetType,
	BaseType,
	TargetParamsType> CorrectionVectorSkeletonType;
	typedef typename BasisType::QnType QnType;

	TargetingCorrectionVector(const LeftRightSuperType& lrs,
	                          const ModelType& model,
	                          const WaveFunctionTransfType& wft,
	                          const QnType&,
	                          InputValidatorType& ioIn)
	    : BaseType(lrs,model,wft,1),
	      tstStruct_(ioIn, model),
	      ioIn_(ioIn),
	      progress_("TargetingCorrectionVector"),
	      gsWeight_(1.0),
	      correctionEnabled_(false),
	      paramsForSolver_(ioIn,"DynamicDmrg"),
	      skeleton_(ioIn_, tstStruct_, model, lrs, this->common().aoe().energy())
	{
		if (!wft.isEnabled())
			err("TargetingCorrectionVector needs wft\n");
	}

	SizeType sites() const { return tstStruct_.sites(); }

	SizeType targets() const { return 4; }

	RealType weight(SizeType i) const
	{
		assert(i < weight_.size());
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (!correctionEnabled_) return 1.0;
		return gsWeight_;
	}

	SizeType size() const
	{
		if (!correctionEnabled_) return 0;
		return BaseType::size();
	}

	void evolve(RealType Eg,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType& block2,
	            SizeType loopNumber)
	{
		if (block1.size()!=1 || block2.size()!=1) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "evolve only blocks of one site supported\n";
			throw PsimagLite::RuntimeError(str.c_str());
		}

		SizeType site = block1[0];
		evolve(Eg,direction,site,loopNumber);

		this->common().printNormsAndWeights(gsWeight_, weight_);

		//corner case
		SizeType numberOfSites = this->lrs().super().block().size();
		SizeType site2 = numberOfSites;

		if (site == 1 && direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
			site2 = 0;
		if (site == numberOfSites - 2 &&
		        direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			site2 = numberOfSites - 1;
		if (site2 == numberOfSites) return;
		evolve(Eg, direction, site2, loopNumber);
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, block, prefix);
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().readGSandNGSTs(io, prefix);
		setWeights();
	}

private:

	void evolve(RealType Eg,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		VectorWithOffsetType phiNew;
		SizeType count = this->common().aoe().getPhi(&phiNew,
		                                             Eg,
		                                             direction,
		                                             site,
		                                             loopNumber,
		                                             tstStruct_);

		if (direction != ProgramGlobals::DirectionEnum::INFINITE) {
			correctionEnabled_=true;
			typename PsimagLite::Vector<SizeType>::Type block1(1,site);
			addCorrection(direction,block1);
		}

		if (count==0) return;

		this->common().aoe().targetVectors(1) = phiNew;
		skeleton_.calcDynVectors(this->common().aoe().targetVectors(1),
		                         this->common().aoe().targetVectors(2),
		                         this->common().aoe().targetVectors(3));

		setWeights();

		VectorSizeType block(1, site);
		bool doBorderIfBorder = false;
		this->common().cocoon(block, direction, doBorderIfBorder);
	}

	void setWeights()
	{
		gsWeight_ = tstStruct_.gsWeight();

		RealType sum  = 0;
		weight_.resize(this->common().aoe().targetVectors().size());
		for (SizeType r=1;r<weight_.size();r++) {
			weight_[r] = 1;
			sum += weight_[r];
		}

		for (SizeType r=0;r<weight_.size();r++) weight_[r] *= (1.0 - gsWeight_)/sum;
	}

	RealType dynWeightOf(VectorType& v,const VectorType& w) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<v.size();i++) {
			RealType tmp = PsimagLite::real(v[i]*w[i]);
			sum += tmp*tmp;
		}
		return sum;
	}

	void addCorrection(ProgramGlobals::DirectionEnum direction,const BlockType& block1)
	{
		if (tstStruct_.correctionA() == 0) return;
		weight_.resize(1);
		weight_[0]=tstStruct_.correctionA();
		this->common().computeCorrection(direction,block1);
		gsWeight_ = 1.0-weight_[0];
	}

	TargetParamsType tstStruct_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	bool correctionEnabled_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	CorrectionVectorSkeletonType skeleton_;
}; // class TargetingCorrectionVector
} // namespace
/*@}*/
#endif // TARGETING_CORRECTION_VECTOR_H

