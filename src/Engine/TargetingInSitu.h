/*
Copyright (c) 2009, 2017, UT-Battelle, LLC
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

/*! \file TargetingInSitu.h
 *
 * Implements the multi-targeting for measuring
 * 2-point correlations in situ
 *
 */

#ifndef DMRG_TARGETING_IN_SITU_H
#define DMRG_TARGETING_IN_SITU_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ParametersForSolver.h"
#include "TargetParamsCorrelations.h"
#include "VectorWithOffsets.h"
#include "TargetingBase.h"
#include <cassert>
#include "Concurrency.h"
#include "Parallelizer.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingInSitu : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef TargetParamsCorrelations<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename BaseType::TargetingCommonType::VectorVectorWithOffsetType
	VectorVectorWithOffsetType;
	typedef typename BaseType::TargetingCommonType::BraketType BraketType;

	enum StageEnum {DISABLED,CONVERGING};

	static SizeType const PRODUCT = TargetParamsType::PRODUCT;
	static SizeType const SUM = TargetParamsType::SUM;

	TargetingInSitu(const LeftRightSuperType& lrs,
	                const ModelType& model,
	                const WaveFunctionTransfType& wft,
	                const SizeType&,
	                InputValidatorType& io)
	    : BaseType(lrs,model,wft,0),
	      tstStruct_(io,model),
	      wft_(wft),
	      progress_("TargetingInSitu"),
	      gsWeight_(tstStruct_.gsWeight()),
	      paramsForSolver_(io,"CorrelationsDmrg"),
	      weightForContinuedFraction_(0),
	      stage_(DISABLED)
	{
		if (tstStruct_.sites() == 0)
			err("TargetingInSitu needs at least one TSPSite\n");

		if (!wft.isEnabled())
			err("TargetingInSitu needs an enabled wft\n");

		this->common().init(&tstStruct_,tstStruct_.sites());

		setWeights();
	}

	RealType weight(SizeType i) const
	{
		return (this->common().targetVectors(i).size() == 0) ?
		            0 : weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeight_;
	}

	SizeType size() const
	{
		return (stage_ == CONVERGING) ?
		            this->common().targetVectors().size() : 0;
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
		SizeType numberOfSites = this->lrs().super().block().size();
		if (site>1 && site<numberOfSites-2) return;
		// //corner case
		//SizeType x = (site==1) ? 0 : numberOfSites-1;
		//evolve(Eg,direction,x,loopNumber);
	}

	void write(const typename PsimagLite::Vector<SizeType>::Type& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix,
	           SizeType counter) const
	{
		this->common().write(io, block, prefix, counter);

		assert(block.size()==1);

		SizeType type = tstStruct_.type();
		int fermionSign = this->common().findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		if (ab_.size()<2) return;
		typename PostProcType::ParametersType params = paramsForSolver_;
		params.Eg = this->common().energy();
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		if (tstStruct_.aOperators()[0].fermionSign>0) s2 *= s;

		PostProcType cf(ab_, params);
		this->common().writeNGSTs(block, io, cf);
	}

	void read(const PsimagLite::String& f)
	{
		this->common().template read<TimeSerializerType>(f);
	}

private:

	void evolve(RealType,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		if (direction == ProgramGlobals::INFINITE) return;

		// see if operator at site has been applied and result put into targetVectors[site]
		// if no apply operator at site and add into targetVectors[site]
		// also wft everything
		this->common().wftAll(site);

		int indexOfOperator = 0;
		SizeType start = 0;
		while ((indexOfOperator = tstStruct_.findIndexOfSite(site, start)) >= 0) {
			this->common().applyOneOperator(loopNumber,
			                                indexOfOperator,
			                                site,
			                                this->common().targetVectors(indexOfOperator),
			                                this->common().psi(),
			                                direction);
			start = indexOfOperator + 1;
		}

		typename PsimagLite::Vector<SizeType>::Type block(1,site);
		cocoon(block,direction);

		if (stage_ == CONVERGING) return;

		for (SizeType i = 0; i < this->common().targetVectors().size(); ++i) {
			if (this->common().targetVectors(i).size() != 0) {
				stage_ = CONVERGING;
				break;
			}
		}
	}

	void setWeights()
	{
		RealType sum  = 0;
		weight_.resize(this->common().targetVectors().size());
		for (SizeType r=0;r<weight_.size();r++) {
			weight_[r] = 1.0;
			sum += weight_[r];
		}

		for (SizeType r=0;r<weight_.size();r++) weight_[r] *=(1.0-gsWeight_)/sum;
	}

	void cocoon(const BlockType& block,
	            ProgramGlobals::DirectionEnum direction) const
	{
		assert(block.size() > 0);
		std::cout<<"TargetingInSitu\n";
		std::cout<<"-------------&*&*&* In-situ measurements start\n";
		typename BraketType::VectorStringType vecStr;
		PsimagLite::split(vecStr, this->model().params().insitu, ",");

		for (SizeType i = 0; i < vecStr.size(); ++i) {
			BraketType braket(this->model(), vecStr[i]);
			this->common().calcBracket(direction, block[0], braket);
		}

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	TargetParamsType tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	TridiagonalMatrixType ab_;
	RealType weightForContinuedFraction_;
	StageEnum stage_;
}; // class TargetingInSitu
} // namespace
/*@}*/
#endif // DMRG_TARGETING_IN_SITU_H

