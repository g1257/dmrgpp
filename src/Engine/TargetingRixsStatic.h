/*
Copyright (c) 2009-2016, UT-Battelle, LLC
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

/*! \file TargetingRixsStatic.h
 *
 * Implements the targeting required by
 * RIXS Static
 * tv[3*site] = A^\dagger_{site}|gs>
 * tv[3*site+1] = Imaginary of (w*-Htilde+i\eta)^{-1}A^\dagger_{site}|gs>
 * tv[3*site+2] = Real of (w*-Htilde+i\eta)^{-1}A^\dagger_{site}|gs>
 *
 */

#ifndef TARGETING_RIXS_STATIC_H
#define TARGETING_RIXS_STATIC_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TargetParamsCorrectionVector.h"
#include "VectorWithOffsets.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "TimeSerializer.h"
#include "FreqEnum.h"
#include "CorrectionVectorSkeleton.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingRixsStatic : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;

public:

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
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,TargetParamsType>
	CorrectionVectorFunctionType;
	typedef ParallelTriDiag<ModelType,LanczosSolverType,VectorWithOffsetType>
	ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename BaseType::InputSimpleOutType InputSimpleOutType;
	typedef CorrectionVectorSkeleton<LanczosSolverType,
	VectorWithOffsetType,
    BaseType,
    TargetParamsType> CorrectionVectorSkeletonType;

	enum StageEnum {STAGE_DISABLED, STAGE_OPERATOR, STAGE_STATIC1, STAGE_STATIC2};

	static SizeType const PRODUCT = TargetParamsType::PRODUCT;
	static SizeType const SUM = TargetParamsType::SUM;

	TargetingRixsStatic(const LeftRightSuperType& lrs,
	                    const ModelType& model,
	                    const WaveFunctionTransfType& wft,
	                    const SizeType&,
	                    InputValidatorType& ioIn)
	    : BaseType(lrs,model,wft,1),
	      stage_(STAGE_DISABLED),
	      tstStruct_(ioIn,model),
	      ioIn_(ioIn),
	      progress_("TargetingRixsStatic"),
	      gsWeight_(1.0),
	      paramsForSolver_(ioIn,"DynamicDmrg"),
	      skeleton_(ioIn_,tstStruct_,model,lrs,this->common().energy())
	{
		SizeType numberOfSites = model.geometry().numberOfSites();
		this->common().init(&tstStruct_,3*numberOfSites);
		if (!wft.isEnabled())
			throw PsimagLite::RuntimeError("TargetingRixsStatic needs wft\n");
	}

	RealType weight(SizeType i) const
	{
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeight_;
	}

	SizeType size() const
	{
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
		SizeType numberOfSites = this->lrs().super().block().size();
		if (site>1 && site<numberOfSites-2) return;
		if (site == 1 && direction == ProgramGlobals::EXPAND_SYSTEM) return;
		//corner case
		//		SizeType x = (site==1) ? 0 : numberOfSites-1;
		//		evolve(Eg,direction,x,loopNumber);

		// skeleton_.printNormsAndWeights();
	}

	void print(InputSimpleOutType& ioOut) const
	{
		ioOut.print("TARGETSTRUCT",tstStruct_);
		PsimagLite::OstringStream msg;
		msg<<"PSI\n";
		msg<<(*this);
		ioOut.print(msg.str());
	}

	void save(const VectorSizeType& block,
	          PsimagLite::IoSimple::Out& io) const
	{
		assert(block.size() > 0);

		SizeType marker = (this->common().noStageIs(STAGE_DISABLED)) ? 1 : 0;
		TimeSerializerType ts(this->common().currentTime(),
		                      block[0],
		                      this->common().targetVectors(),
		                      marker);
		ts.save(io);
		const VectorWithOffsetType& v = this->common().targetVectors()[3];
		this->common().psi().save(io,"PSI");

		for (SizeType i = 0; i < v.sectors(); ++i) {
			SizeType m = v.sector(i);
			SizeType state = this->lrs().super().partition(m);
			SizeType qn = this->lrs().super().pseudoEffectiveNumber(state);
			PsimagLite::String str("PSI.QN." + ttos(m) + "=" + ttos(qn));
			io.printline(str);
		}
	}

	void load(const PsimagLite::String& f)
	{
		this->common().template load<TimeSerializerType>(f);
	}

private:

	// tv[i] = A_site |gs> if i%3 == 0. site = i/3.
	// tv[i+1] = imaginary cv for tv[i]
	// tv[i+2] = real      cv for tv[i]
	void evolve(RealType Eg,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		if (direction == ProgramGlobals::INFINITE) return;
		SizeType numberOfSites = this->lrs().super().block().size();
		SizeType indexOfOperator = 0;

		// see if operator at site has been applied and result put into targetVectors[site]
		// if no apply operator at site and add into targetVectors[site]
		// also wft everything
		if (stage_ != STAGE_STATIC2) {
			VectorSizeType indexForOperators(this->common().targetVectors().size(), 0);
			this->common().wftAll(indexForOperators, site, direction);
			this->common().applyOneOperator(loopNumber,
			                                indexOfOperator,
			                                site,
			                                this->common().targetVectors(3*site),
			                                direction);

		} else {
			for (SizeType s = 0; s < numberOfSites; ++s) {
				const VectorWithOffsetType& src = this->common().targetVectors(3*s);
				if (src.size() == 0) continue;
				VectorWithOffsetType phiNew;
				if (tstStruct_.useQns()) phiNew.populateFromQns(this->common().nonZeroQns(),
				                                                this->lrs().super());
				this->common().wftOneVector(phiNew,
				                            src,
				                            indexOfOperator,
				                            site,
				                            direction);
				this->common().targetVectors(3*s) = phiNew;
			}

			doCorrectionVector(direction, site);
		}

		// typename PsimagLite::Vector<SizeType>::Type block(1,site);
		// skeleton_.cocoon(block,direction);

		if (stage_ == STAGE_STATIC2) return;

		SizeType counter = 0;
		for (SizeType i = 0; i < this->common().targetVectors().size(); ++i) {
			if (i%3 != 0) continue;
			if (this->common().targetVectors(i).size() != 0) {
				if (stage_ != STAGE_STATIC2) stage_ = STAGE_STATIC1;
				counter++;
			}
		}

		// -2 here because borders not done (needs FIXING FIXME TODO)
		if (counter == numberOfSites-2 && stage_ == STAGE_STATIC1)
			stage_ = STAGE_STATIC2;
	}

	void doCorrectionVector(ProgramGlobals::DirectionEnum direction,
	                        SizeType site)
	{
		assert(stage_ == STAGE_STATIC2);
		SizeType numberOfSites = this->lrs().super().block().size();
		for (SizeType s = 0; s < numberOfSites; ++s)
			skeleton_.calcDynVectors(this->common().targetVectors(3*s),
			                         this->common().targetVectors(3*s+1),
			                         this->common().targetVectors(3*s+2),
			                         direction,
			                         site);
		setWeights();
	}

	void setWeights()
	{
		gsWeight_ = tstStruct_.gsWeight();

		RealType sum  = 0;
		weight_.resize(this->common().targetVectors().size());
		for (SizeType r=1;r<weight_.size();r++) {
			weight_[r] = 1;
			sum += weight_[r];
		}

		for (SizeType r=0;r<weight_.size();r++) weight_[r] *= (1.0 - gsWeight_)/sum;
	}

	void printNormsAndWeights() const
	{
		if (this->common().allStages(STAGE_DISABLED)) return;

		PsimagLite::OstringStream msg;
		msg<<"gsWeight="<<gsWeight_<<" weights= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg<<weight_[i]<<" ";
		progress_.printline(msg,std::cout);

		PsimagLite::OstringStream msg2;
		msg2<<"gsNorm="<<norm(this->common().psi())<<" norms= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg2<<this->common().normSquared(i)<<" ";
		progress_.printline(msg2,std::cout);
	}

	StageEnum stage_;
	TargetParamsType tstStruct_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	CorrectionVectorSkeletonType skeleton_;
}; // class TargetingRixsStatic

template<typename LanczosSolverType, typename VectorWithOffsetType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingRixsStatic<LanczosSolverType,VectorWithOffsetType>&)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // TARGETING_RIXS_STATIC_H
