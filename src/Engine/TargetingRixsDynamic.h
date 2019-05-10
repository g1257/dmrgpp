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
#include "TargetParamsTimeStep.h"
#include "VectorWithOffsets.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "FreqEnum.h"
#include "CorrectionVectorSkeleton.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingRixsDynamic : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

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
	typedef TargetParamsTimeStep<ModelType> TargetParams2Type;
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
	typedef typename BasisType::QnType QnType;
	typedef typename TargetParamsType::BaseType::AlgorithmEnum AlgorithmEnumType;
	typedef typename TargetingCommonType::StageEnumType StageEnumType;

	TargetingRixsDynamic(const LeftRightSuperType& lrs,
	                     const ModelType& model,
	                     const WaveFunctionTransfType& wft,
	                     const QnType&,
	                     InputValidatorType& ioIn)
	    : BaseType(lrs,model,wft,1),
	      tstStruct_(ioIn, "TargetingRixsDynamic", model),
	      tstStruct2_(nullptr),
	      ioIn_(ioIn),
	      progress_("TargetingRixsDynamic"),
	      gsWeight_(1.0),
	      paramsForSolver_(ioIn,"DynamicDmrg"),
	      skeleton_(ioIn_,tstStruct_,model,lrs,this->common().aoe().energy()),
	      applied_(false),
	      appliedFirst_(false),
	      usesCheby_(tstStruct_.algorithm() == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV)
	{
		if (!wft.isEnabled())
			err("TargetingRixsDynamic needs wft\n");

		if (!usesCheby_) {
			return; // early exit here
		}

		tstStruct2_ = new TargetParams2Type(ioIn, "TargetingRixsDynamic", model);
		times_.resize(tstStruct2_->timeSteps());

		RealType tau = tstStruct2_->tau();
		SizeType n = times_.size();
		for (SizeType i = 0; i < n; ++i)
			times_[i] = i*tau/(n - 1);

		this->common().aoe().initTimeVectors(*tstStruct2_, times_, ioIn);
	}

	~TargetingRixsDynamic()
	{
		delete tstStruct2_;
		tstStruct2_ = nullptr;
	}

	SizeType sites() const { return (usesCheby_) ? tstStruct2_->sites() : tstStruct_.sites(); }

	SizeType targets() const { return 12; }

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
		if (!applied_ && appliedFirst_) return 8;

		const AlgorithmEnumType algo = tstStruct_.algorithm();
		SizeType tenOrTwelve = (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV) ? 12
		                                                                                      : 10;
		return (applied_) ? tenOrTwelve : 6;
	}

	// tv[6] = A^\dagger_{site} |tv[1]>
	// tv[7] = A^\dagger_{site} |tv[2]>

	// tv[8] = Re of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[6]>
	// - Im of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[7]>
	// tv[9] = Im of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[6]>
	// + Re of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[7]>
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

		if (direction == ProgramGlobals::DirectionEnum::INFINITE) return;

		SizeType max = tstStruct_.sites();

		if (max > 2)
			err("You cannot apply more than 2 operators (only SUM is allowed)\n");

		SizeType site = block1[0];

		this->common().aoe().wftSome(site, 0, 8);

		const AlgorithmEnumType algo = tstStruct_.algorithm();
		if (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV) {
			// just to set the stage and currenttime
			this->common().aoe().getPhi(0, Eg, direction, site, loopNumber, *tstStruct2_);
		} else if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			this->common().aoe().wftSome(site, 8, this->common().aoe().targetVectors().size());
		} else {
			assert(false);
		}

		if (!applied_) {
			if (max == 1)
				doMax1(site, direction, loopNumber);

			if (max == 2)
				doMax2(site, direction, loopNumber);
		}

		calcDynVectors(Eg, direction, block1); // WEIGHTS ARE SET IN calcDynVectors()

		cocoon(site, direction);
	}

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, block, prefix);
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().read(io, prefix);

		TimeSerializerType ts(io, prefix);
		SizeType n = ts.numberOfVectors();
		if (n != 6)
			err("TargetingRixsDynamic: number of TVs must be 6\n");

		for (SizeType site = 0; site < 6; ++site)
			this->common().aoe().targetVectors(site) = ts.vector(site);
	}

private:

	void doMax1(SizeType site,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType loopNumber)
	{
		if (site == tstStruct_.sites(0)) {

			ComplexOrRealType densCim = this->common().rixsCocoon(direction,site,1,0,false);
			std::cout<<site<<" "<<densCim<<" 0"; // 0 here is the currentTime
			std::cout<<" <P1|P0> 1\n";   // 1 here is the "superdensity"

			ComplexOrRealType densCre = this->common().rixsCocoon(direction,site,2,0,false);
			std::cout<<site<<" "<<densCre<<" 0"; // 0 here is the currentTime
			std::cout<<" <P2|P0> 1\n";   // 1 here is the "superdensity"

			ComplexOrRealType densjim = this->common().rixsCocoon(direction,site,4,3,false);
			std::cout<<site<<" "<<densjim<<" 0"; // 0 here is the currentTime
			std::cout<<" <P4|P3> 1\n";   // 1 here is the "superdensity"

			ComplexOrRealType densjre = this->common().rixsCocoon(direction,site,5,3,false);
			std::cout<<site<<" "<<densjre<<" 0"; // 0 here is the currentTime
			std::cout<<" <P5|P3> 1\n";   // 1 here is the "superdensity"

			VectorWithOffsetType tmpV1;
			SizeType indexOfOperator = 0;

			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV1, // phiNew
			           this->common().aoe().targetVectors(1), // src1 apply op on Im|alpha(C)>
			           direction);

			if (tmpV1.size() > 0)
				addFactor(tmpV1, this->common().aoe().psi(), densCim);

			if (tmpV1.size() > 0)
				this->common().aoe().targetVectors(6) = tmpV1;

			VectorWithOffsetType tmpV2;

			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->common().aoe().targetVectors(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0)
				addFactor(tmpV2, this->common().aoe().psi(), densCre);

			if (tmpV2.size() > 0) {
				this->common().aoe().targetVectors(7) = tmpV2;
				applied_ = true;
				PsimagLite::OstringStream msg;
				msg<<"Applied";
				progress_.printline(msg, std::cout);
			}
		}
	}

	void doMax2(SizeType site,
	            ProgramGlobals::DirectionEnum direction,
	            SizeType loopNumber)
	{

		if (site == tstStruct_.sites(0)) {
			VectorWithOffsetType tmpV1;
			SizeType indexOfOperator = 0;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV1, // phiNew
			           this->common().aoe().targetVectors(1), // src1 apply op on Im|alpha(C)>
			           direction);

			if (tmpV1.size() > 0)
				this->common().aoe().targetVectors(6) = tmpV1;

			VectorWithOffsetType tmpV2;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->common().aoe().targetVectors(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0) {
				this->common().aoe().targetVectors(7) = tmpV2;
				applied_ = false;
				appliedFirst_ = true;
				PsimagLite::OstringStream msg;
				msg<<"First Operator Applied";
				progress_.printline(msg, std::cout);
			}
		}

		if (site == tstStruct_.sites(1)) {

			ComplexOrRealType densCim = this->common().rixsCocoon(direction,site,1,0,false);
			std::cout<<site<<" "<<densCim<<" 0"; // 0 here is the currentTime
			std::cout<<" <P1|P0> 1\n";   // 1 here is the "superdensity"

			ComplexOrRealType densCre = this->common().rixsCocoon(direction,site,2,0,false);
			std::cout<<site<<" "<<densCre<<" 0"; // 0 here is the currentTime
			std::cout<<" <P2|P0> 1\n";   // 1 here is the "superdensity"

			ComplexOrRealType densjim = this->common().rixsCocoon(direction,site,4,3,false);
			std::cout<<site<<" "<<densjim<<" 0"; // 0 here is the currentTime
			std::cout<<" <P4|P3> 1\n";   // 1 here is the "superdensity"

			ComplexOrRealType densjre = this->common().rixsCocoon(direction,site,5,3,false);
			std::cout<<site<<" "<<densjre<<" 0"; // 0 here is the currentTime
			std::cout<<" <P5|P3> 1\n";   // 1 here is the "superdensity"

			VectorWithOffsetType tmpV1;
			SizeType indexOfOperator = 1;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV1, // phiNew
			           this->common().aoe().targetVectors(1), // src1 apply op on Im|alpha(C)>
			           direction);

			if (tmpV1.size() > 0)
				addFactor(tmpV1, this->common().aoe().psi(), densCim);

			if (tmpV1.size() > 0)
				this->common().aoe().targetVectors(6) += tmpV1;

			VectorWithOffsetType tmpV2;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->common().aoe().targetVectors(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0)
				addFactor(tmpV2, this->common().aoe().psi(), densCre);

			if (tmpV2.size() > 0) {
				this->common().aoe().targetVectors(7) += tmpV2;
				applied_ = true;
				PsimagLite::OstringStream msg;
				msg<<"Applied";
				progress_.printline(msg, std::cout);
			}
		}
	}

	void cocoon(SizeType site, ProgramGlobals::DirectionEnum direction)
	const
	{
		if (!usesCheby_) {
			ComplexOrRealType rr =
			        this->common().rixsCocoon(direction,site,9,5,true);
			ComplexOrRealType ri =
			        this->common().rixsCocoon(direction,site,9,4,true);
			ComplexOrRealType ir =
			        this->common().rixsCocoon(direction,site,8,5,true);
			ComplexOrRealType ii =
			        this->common().rixsCocoon(direction,site,8,4,true);

			const RealType time = this->common().aoe().time();
			std::cout<<site<<" "<<(ri-ir)<<" "<<time; // time here is the currentTime
			std::cout<<" <gs|A|P2> 1\n";   // 1 here is the "superdensity"
			std::cout<<site<<" "<<(rr+ii)<<" "<<time; // time here is the currentTime
			std::cout<<" <gs|A|P3> 1\n";   // 1 here is the "superdensity"
			return;
		} else {
			ComplexOrRealType rr =
			        this->common().rixsCocoon(direction,site,10,5,true);
			ComplexOrRealType ri =
			        this->common().rixsCocoon(direction,site,10,4,true);
			ComplexOrRealType ir =
			        this->common().rixsCocoon(direction,site,8,5,true);
			ComplexOrRealType ii =
			        this->common().rixsCocoon(direction,site,8,4,true);

			const RealType time = this->common().aoe().time();
			std::cout<<site<<" "<<(ri-ir)<<" "<<time; // time here is the currentTime
			std::cout<<" <gs|A|P2> 1\n";   // 1 here is the "superdensity"
			std::cout<<site<<" "<<(rr+ii)<<" "<<time; // time here is the currentTime
			std::cout<<" <gs|A|P3> 1\n";   // 1 here is the "superdensity"

			for (SizeType i = 6; i < 12; ++i)
				std::cout<<"norm2("<<i<<")= "<<this->common().normSquared(i)<<"\n";
		}
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

	void calcDynVectors(RealType Eg,
	                    ProgramGlobals::DirectionEnum direction,
	                    const VectorSizeType& block1)
	{

		if (!applied_ && appliedFirst_) {
			setWeights(8);
			return;
		}

		if (!applied_) {
			setWeights(6);
			return;
		}

		const AlgorithmEnumType algo = tstStruct_.algorithm();

		if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			skeleton_.calcDynVectors(this->common().aoe().targetVectors(6),
			                         this->common().aoe().targetVectors(7),
			                         this->common().aoe().targetVectors(8),
			                         this->common().aoe().targetVectors(9));
			setWeights(10);
		} else if (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV) {
			VectorSizeType indices{6, 8, 9};
			calcChebyVectors(indices, Eg, direction, block1);
			VectorSizeType indices2{7, 10, 11};
			calcChebyVectors(indices2, Eg, direction, block1);
			setWeights(12);
		} else {
			assert(false);
		}
	}

	void calcChebyVectors(const VectorSizeType& indices,
	                      RealType Eg,
	                      ProgramGlobals::DirectionEnum direction,
	                      const VectorSizeType& block1)
	{
		assert(indices.size() >= 3);
		bool allOperatorsApplied = (this->common().aoe().noStageIs(StageEnumType::DISABLED) &&
		                            this->common().aoe().noStageIs(StageEnumType::OPERATOR));

		const VectorWithOffsetType& v0 = this->common().aoe().targetVectors(indices[0]);
		this->common().chebyshev(indices,
		                         Eg,
		                         v0,
		                         direction,
		                         allOperatorsApplied,
		                         block1,
		                         *tstStruct2_);
	}

	void applyOneOp(SizeType loopNumber,
	                SizeType indexOfOperator,
	                SizeType site,
	                VectorWithOffsetType& dest,
	                const VectorWithOffsetType& src,
	                ProgramGlobals::DirectionEnum direction)
	{
		if (usesCheby_)
			this->common().aoe().applyOneOperator(loopNumber,
		                                          indexOfOperator,
		                                          site,
		                                          dest, // phiNew
		                                          src, // src1
		                                          direction,
		                                          *tstStruct2_);
		else
			this->common().aoe().applyOneOperator(loopNumber,
		                                          indexOfOperator,
		                                          site,
		                                          dest, // phiNew
		                                          src, // src1
		                                          direction,
		                                          tstStruct_);
	}

	void setWeights(SizeType n)
	{
		gsWeight_ = tstStruct_.gsWeight();

		RealType sum  = n;
		weight_.resize(n, 1);

		for (SizeType r=0;r<weight_.size();r++) weight_[r] = (1.0 - gsWeight_)/sum;
	}

	TargetParamsType tstStruct_;
	TargetParams2Type* tstStruct2_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	VectorRealType times_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	CorrectionVectorSkeletonType skeleton_;
	bool applied_;
	bool appliedFirst_;
	bool usesCheby_;
}; // class TargetingRixsDynamic
} // namespace
/*@}*/
#endif // TARGETING_RIXS_DYNAMIC_H
