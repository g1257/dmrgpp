/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
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
 * We read from static tv[i] --> tv[i]
 *
 * the correction vectors are imag  --> tv[8]
 *                            real  --> tv[9]
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
	      skeleton_(ioIn_,tstStruct_,model,lrs,this->common().energy()),
	      applied_(false)
	{
		//		SizeType numberOfSites = model.geometry().numberOfSites();
		this->common().init(&tstStruct_,12);
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
		return (applied_) ? 10 : 6;
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

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix,
	           SizeType counter) const
	{
		this->common().write(io, block, prefix, counter);
		this->common().writeNGSTs(block, io);
	}

	void read(const PsimagLite::String& f)
	{
		PsimagLite::IoSelector::In io(f);

		TimeSerializerType ts(io, PsimagLite::IoSimple::In::LAST_INSTANCE);
		SizeType n = ts.numberOfVectors();
		if (n != 6)
			err("TargetingRixsDynamic: number of TVs must be 6\n");

		for (SizeType site = 0; site < 6; ++site) {
			this->common().targetVectors(site) = ts.vector(site);
		}

		this->common().template read<TimeSerializerType>(f,0);
	}

private:

	// tv[6] = A^\dagger_{site} |tv[1]>
	// tv[7] = A^\dagger_{site} |tv[2]>

	// tv[8] = Re of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[6]>
	// - Im of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[7]>
	// tv[9] = Im of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[6]>
	// + Re of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[7]>

	void evolve(RealType,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType site,
	            SizeType loopNumber)
	{
		if (direction == ProgramGlobals::INFINITE) return;

		SizeType indexOfOperator = 0;
		this->common().wftAll(site);

		ComplexOrRealType densCre = this->common().rixsCocoon(direction,site,1,0,false);
		std::cout<<site<<" "<<densCre<<" 0"; // 0 here is the currentTime
		std::cout<<" <P1|P0> 1\n";   // 1 here is the "superdensity"

		ComplexOrRealType densCim = this->common().rixsCocoon(direction,site,2,0,false);
		std::cout<<site<<" "<<densCim<<" 0"; // 0 here is the currentTime
		std::cout<<" <P2|P0> 1\n";   // 1 here is the "superdensity"

		ComplexOrRealType densjre = this->common().rixsCocoon(direction,site,4,3,false);
		std::cout<<site<<" "<<densjre<<" 0"; // 0 here is the currentTime
		std::cout<<" <P4|P3> 1\n";   // 1 here is the "superdensity"

		ComplexOrRealType densjim = this->common().rixsCocoon(direction,site,5,3,false);
		std::cout<<site<<" "<<densjim<<" 0"; // 0 here is the currentTime
		std::cout<<" <P5|P3> 1\n";   // 1 here is the "superdensity"

		if (!applied_) {
			if (site == tstStruct_.sites(0)) {

				VectorWithOffsetType tmpV1;
				this->common().applyOneOperator(loopNumber,
				                                indexOfOperator,
				                                site,
				                                tmpV1, // phiNew
				                                this->common().targetVectors(1), // src1
				                                direction);

				addFactor(tmpV1, this->common().psi(), densCre);

				if (tmpV1.size() > 0)
					this->common().targetVectors(6) = tmpV1;

				VectorWithOffsetType tmpV2;
				this->common().applyOneOperator(loopNumber,
				                                indexOfOperator,
				                                site,
				                                tmpV2,                            // phiNew
				                                this->common().targetVectors(2), // src1
				                                direction);

				addFactor(tmpV2, this->common().psi(), densCim);

				if (tmpV2.size() > 0) {
					this->common().targetVectors(7) = tmpV2;
					applied_ = true;
					PsimagLite::OstringStream msg;
					msg<<"Applied";
					progress_.printline(msg, std::cout);
				}
			}
		}
		calcDynVectors();

		// WEIGHTS ARE SET IN calcDynVectors()

		ComplexOrRealType rr = this->common().rixsCocoon(direction,site,8,4,true);
		ComplexOrRealType ri = this->common().rixsCocoon(direction,site,8,5,true);
		ComplexOrRealType ir = this->common().rixsCocoon(direction,site,9,4,true);
		ComplexOrRealType ii = this->common().rixsCocoon(direction,site,9,5,true);

		std::cout<<site<<" "<<(ri-ir)<<" 0"; // 0 here is the currentTime
		std::cout<<" <gs|A|P2> 1\n";   // 1 here is the "superdensity"
		std::cout<<site<<" "<<(rr+ii)<<" 0"; // 0 here is the currentTime
		std::cout<<" <gs|A|P3> 1\n";   // 1 here is the "superdensity"
	}

	void addFactor(VectorWithOffsetType& phiNew,
	               const VectorWithOffsetType& psiSrc2,
	               ComplexOrRealType factor) const
	{
		// CHECK if psiSrc2 and phiNew have the same offset!
		if (psiSrc2.offset(0) == phiNew.offset(0))
			phiNew += (-factor)*psiSrc2;

		RealType norma = norm(phiNew);
		if (norma<1e-6) {
			PsimagLite::OstringStream msg2;
			msg2<<"Norm of phi is zero\n";
			progress_.printline(msg2,std::cout);
		}
	}

	void calcDynVectors()
	{
		if (!applied_) {
			setWeights(6);
			return;
		}
		skeleton_.calcDynVectors(this->common().targetVectors(6),
		                         this->common().targetVectors(7),
		                         this->common().targetVectors(8),
		                         this->common().targetVectors(9));
		setWeights(10);
	}

	void setWeights(SizeType n)
	{
		gsWeight_ = tstStruct_.gsWeight();

		RealType sum  = n;
		weight_.resize(n, 1);

		for (SizeType r=0;r<weight_.size();r++) weight_[r] = (1.0 - gsWeight_)/sum;
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
	bool applied_;
}; // class TargetingRixsDynamic
} // namespace
/*@}*/
#endif // TARGETING_RIXS_DYNAMIC_H
