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

/*! \file TargetingRixsDynamic.h
 *
 * Implements the targeting required by
 * RIXS Dynamic
 *
 * Must be restarted from RIXS Static
 *
 * We read from static tv[3*site+1] --> tv[2*site]
 *                     tv[3*site+2] --> tv[2*site+1]
 * the correction vectors are imag  --> tv[2*N]
 *                            real  --> tv[2*N+1]
 *
 */

#ifndef TARGETING_RIXS_DYNAMIC_H
#define TARGETING_RIXS_DYNAMIC_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TargetParamsCorrectionVector.h"
#include "VectorWithOffsets.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "FreqEnum.h"
#include "CorrectionVectorSkeleton.h"
#include "TimeSerializer.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingRixsDynamic : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

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
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
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

	enum {DISABLED,OPERATOR,CONVERGING};

	static SizeType const PRODUCT = TargetParamsType::PRODUCT;
	static SizeType const SUM = TargetParamsType::SUM;

	TargetingRixsDynamic(const LeftRightSuperType& lrs,
	                     const ModelType& model,
	                     const WaveFunctionTransfType& wft,
	                     const SizeType&,
	                     InputValidatorType& ioIn)
	    : BaseType(lrs,model,wft,1),
	      tstStruct_(ioIn,model),
	      ioIn_(ioIn),
	      progress_("TargetingRixsDynamic"),
	      gsWeight_(1.0),
	      paramsForSolver_(ioIn,"DynamicDmrg"),
	      skeleton_(ioIn_,tstStruct_,model,lrs,this->common().energy())
	{
		if (tstStruct_.concatenation() != TargetParamsType::DONT_APPLY)
			err("TargetingRixsDynamic needs TSPProductOrSum=dontapply\n");

		SizeType numberOfSites = model.geometry().numberOfSites();
		this->common().init(&tstStruct_,2*numberOfSites+2);
		if (!wft.isEnabled())
			throw PsimagLite::RuntimeError("TargetingRixsDynamic needs wft\n");
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
		skeleton_.save(this->common(),block,io);
	}

	void load(const PsimagLite::String& f)
	{
		typename BaseType::IoInputType io(f);

		TimeSerializerType ts(io,BaseType::IoInputType::LAST_INSTANCE);
		SizeType n = ts.numberOfVectors();
		if (n % 3 != 0)
			err("TargetingRixsDynamic: number of TVs not divisible by 3\n");

		SizeType numberOfSites = n/3;
		for (SizeType site = 0; site < numberOfSites; ++site) {
			this->common().targetVectors(2*site) = ts.vector(3*site + 1);
			this->common().targetVectors(2*site + 1) = ts.vector(3*site + 2);
		}

		this->common().template load<TimeSerializerType>(f,0);
	}

private:

	// tv[2*site] = Imaginary of (w*-Htilde+i\eta)^{-1}A^\dagger_{site}|gs>
	// tv[2*site+1] = Real of (w*-Htilde+i\eta)^{-1}A^\dagger_{site}|gs>
	// tv[2*N] = imaginary cv for (tv[2*center],tv[2*center+1])
	// tv[2*N+1] = real      cv for (tv[2*center],tv[2*center+1])
	void evolve(RealType,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		if (direction == ProgramGlobals::INFINITE) return;

		SizeType indexOfOperator = 0;
		VectorSizeType indexForOperators(this->common().targetVectors().size(), 0);
		SizeType center = tstStruct_.sites(indexOfOperator);
		this->common().wftAll(indexForOperators, site,direction);
		this->common().applyOneOperator(loopNumber,
		                                indexOfOperator,
		                                site,
		                                this->common().targetVectors(2*center),
		                                direction);
		this->common().applyOneOperator(loopNumber,
		                                indexOfOperator,
		                                site,
		                                this->common().targetVectors(2*center+1),
		                                direction);

		calcDynVectors(site,direction);

		SizeType numberOfSites = this->lrs().super().block().size();
		ComplexOrRealType rr = this->common().rixsCocoon(direction,site,2*site,2*numberOfSites);
		ComplexOrRealType ri = this->common().rixsCocoon(direction,site,2*site,2*numberOfSites+1);
		ComplexOrRealType ir = this->common().rixsCocoon(direction,site,2*site+1,2*numberOfSites);
		ComplexOrRealType ii = this->common().rixsCocoon(direction,site,2*site+1,2*numberOfSites+1);

		std::cout<<site<<" "<<(ri+ir)<<" 0"; // 0 here is the currentTime
		std::cout<<" <gs|A|P2> 1\n";   // 1 here is the "superdensity"
		std::cout<<site<<" "<<(rr-ii)<<" 0"; // 0 here is the currentTime
		std::cout<<" <gs|A|P3> 1\n";   // 1 here is the "superdensity"
	}

	void calcDynVectors(SizeType site,
	                    ProgramGlobals::DirectionEnum direction)
	{
		SizeType numberOfSites = this->lrs().super().block().size();
		SizeType center = tstStruct_.sites(0);
		skeleton_.calcDynVectors(this->common().targetVectors(2*center),
		                         this->common().targetVectors(2*center+1),
		                         this->common().targetVectors(2*numberOfSites),
		                         this->common().targetVectors(2*numberOfSites+1),
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
		if (this->common().allStages(DISABLED)) return;

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

	TargetParamsType tstStruct_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	CorrectionVectorSkeletonType skeleton_;
}; // class TargetingRixsDynamic

template<typename LanczosSolverType, typename VectorWithOffsetType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingRixsDynamic<LanczosSolverType,VectorWithOffsetType>&)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // TARGETING_RIXS_DYNAMIC_H
