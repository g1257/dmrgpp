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
 *                     real  --> tv[9]
 *
 */

#ifndef TARGETING_RIXS_DYNAMIC_H
#define TARGETING_RIXS_DYNAMIC_H

#include "ProgressIndicator.h"
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
	typedef typename BaseType::CheckpointType CheckpointType;
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
	                     const CheckpointType& checkPoint,
	                     const WaveFunctionTransfType& wft,
	                     const QnType&,
	                     InputValidatorType& ioIn)
	    : BaseType(lrs,checkPoint,wft,1),
	      tstStruct_(ioIn, "TargetingRixsDynamic", checkPoint.model()),
	      tstStruct2_(nullptr),
	      ioIn_(ioIn),
	      progress_("TargetingRixsDynamic"),
	      gsWeight_(tstStruct_.gsWeight()),
	      paramsForSolver_(ioIn,"DynamicDmrg"),
	      skeleton_(ioIn_, tstStruct_, checkPoint.model(), lrs, this->common().aoe().energy()),
	      applied_(false),
	      appliedFirst_(false)
	{
		firstCall_ = true;

		if (!wft.isEnabled())
			err("TargetingRixsDynamic needs wft\n");

		if (tstStruct_.algorithm() == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			return; // early exit here
		}

		tstStruct2_ = new TargetParams2Type(ioIn, "TargetingRixsDynamic", checkPoint.model());

		RealType tau = tstStruct2_->tau();
		SizeType n = tstStruct2_->times().size();
		if (tstStruct_.algorithm() == TargetParamsType::BaseType::AlgorithmEnum::KRYLOVTIME) {
			if (n != 5)
				err("TargetingRixsDynamic with KrylovTime: number of TimeSteps must be 5\n");
		}

		for (SizeType i = 0; i < n; ++i)
			tstStruct2_->times()[i] = i*tau/(n - 1);

		this->common().aoeNonConst().initTimeVectors(*tstStruct2_, ioIn);
	}

	~TargetingRixsDynamic()
	{
		delete tstStruct2_;
		tstStruct2_ = nullptr;
	}

	SizeType sites() const { return (tstStruct_.algorithm() ==
		                             TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) ?
		            tstStruct_.sites() : tstStruct2_->sites(); }

	SizeType targets() const
	{
		const AlgorithmEnumType algo = tstStruct_.algorithm();
		if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			return 10;
		} else if (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV) {
			return 12;
		} else {
			return 16;
		}
	}

	RealType weight(SizeType i) const
	{
		assert(i < weight_.size());
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeightActual_;
	}

	SizeType size() const
	{
		if (!applied_ && appliedFirst_) return 8;

		SizeType tenOrTwelveOrSixteen;
		const AlgorithmEnumType algo = tstStruct_.algorithm();
		if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			tenOrTwelveOrSixteen = 10;
		} else if (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV) {
			tenOrTwelveOrSixteen = 12;
		} else {
			tenOrTwelveOrSixteen = 16;
		}

		return (applied_) ? tenOrTwelveOrSixteen : 6;
	}

	// tv[6] = A^\dagger_{site} |tv[1]>
	// tv[7] = A^\dagger_{site} |tv[2]>

	// tv[8] = Re of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[6]>
	// - Im of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[7]>
	// tv[9] = Im of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[6]>
	// + Re of (w*-H+i\eta)^{-1} A^\dagger_{site} |tv[7]>
	void evolve(const VectorRealType& energies,
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

		assert(energies.size() > 0);
		const RealType Eg = energies[0];
		SizeType max = tstStruct_.sites();

		if (max > 2)
			err("You cannot apply more than 2 operators (only SUM is allowed)\n");

		SizeType site = block1[0];

		SizeType numberOfSites = this->lrs().super().block().size();
		bool weAreAtBorder = site==0 || site==numberOfSites-1;
		if (!weAreAtBorder)
			this->common().aoeNonConst().wftSome(site, 0, 6);

		const AlgorithmEnumType algo = tstStruct_.algorithm();
		if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			if (!weAreAtBorder)
				this->common().aoeNonConst().wftSome(site, 6, this->common().aoe().tvs());
		} else {
			// just to set the stage and currenttime: CHEBY and KRYLOVTIME
			this->common().aoeNonConst().getPhi(0, Eg, direction, site, loopNumber, *tstStruct2_);
		}

		if (!applied_) {
			if (max == 1)
				doMax1(site, direction, loopNumber);

			if (max == 2 && tstStruct_.concatenation() == TargetParamsType::ConcatEnum::SUM)
				doMax2Sum(site, direction, loopNumber);

			if (max == 2 && tstStruct_.concatenation() == TargetParamsType::ConcatEnum::PRODUCT)
				doMax2Prod(site, direction, loopNumber);
		}

		calcDynVectors(Eg, direction, block1); // WEIGHTS ARE SET IN calcDynVectors()

		cocoon(site, direction);

		this->common().printNormsAndWeights(gsWeight_, weight_);

		//corner case
		SizeType site2 = numberOfSites;

		if (site == 1 && direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
			site2 = 0;
		if (site == numberOfSites - 2 &&
		        direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			site2 = numberOfSites - 1;
		if (site2 == numberOfSites) return;
		BlockType block(1, site2);
		evolve(energies, direction, block, block2, loopNumber);
		setWeights();
	}

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);
		this->common().writeNGSTs(io, prefix, block, "RixsDynamic");
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().read(io, prefix);

		TimeSerializerType ts(io, prefix);
		SizeType n = ts.numberOfVectors();
		if (n != 6)
			err("TargetingRixsDynamic: number of TVs must be 6\n");

		for (SizeType site = 0; site < 6; ++site)
			this->tvNonConst(site) = ts.vector(site);

		setWeights();
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
			           this->tv(1), // src1 apply op on Im|alpha(C)>
			           direction);

			const VectorWithOffsetType& psi00 = this->common().aoe().
			        ensureOnlyOnePsi(__FILE__ + PsimagLite::String("::doMax1"));
			if (tmpV1.size() > 0)
				addFactor(tmpV1, psi00, densCim);

			if (tmpV1.size() > 0)
				this->tvNonConst(6) = tmpV1;

			VectorWithOffsetType tmpV2;

			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phi
			           this->tv(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0)
				addFactor(tmpV2, psi00, densCre);

			if (tmpV2.size() > 0) {
				this->tvNonConst(7) = tmpV2;
				applied_ = true;
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"Applied";
				progress_.printline(msgg, std::cout);
			}
		}
	}

	void doMax2Sum(SizeType site,
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
			           this->tv(1), // src1 apply op on Im|alpha(C)>
			           direction);

			if (tmpV1.size() > 0)
				this->tvNonConst(6) = tmpV1;

			VectorWithOffsetType tmpV2;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->tv(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0) {
				this->tvNonConst(7) = tmpV2;
				applied_ = false;
				appliedFirst_ = true;
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"First Operator Applied";
				progress_.printline(msgg, std::cout);
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
			           this->tv(1), // src1 apply op on Im|alpha(C)>
			           direction);

			const VectorWithOffsetType& psi00 = this->common().aoe().
			        ensureOnlyOnePsi(__FILE__ + PsimagLite::String("::doMax2"));

			if (tmpV1.size() > 0)
				addFactor(tmpV1, psi00, densCim);

			if (tmpV1.size() > 0)
				this->tvNonConst(6) += tmpV1;

			VectorWithOffsetType tmpV2;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->tv(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0)
				addFactor(tmpV2, psi00, densCre);

			if (tmpV2.size() > 0) {
				this->tvNonConst(7) += tmpV2;
				applied_ = true;
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"Applied";
				progress_.printline(msgg, std::cout);
			}
		}
	}

	void doMax2Prod(SizeType site,
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
			           this->tv(1), // src1 apply op on Im|alpha(C)>
			           direction);

			if (tmpV1.size() > 0)
				this->tvNonConst(6) = tmpV1;

			VectorWithOffsetType tmpV2;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->tv(2), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0) {
				this->tvNonConst(7) = tmpV2;
				applied_ = false;
				appliedFirst_ = true;
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"PROD: First Operator Applied";
				progress_.printline(msgg, std::cout);
			}
		}

		if (site == tstStruct_.sites(1)) {

			VectorWithOffsetType tmpV1;
			SizeType indexOfOperator = 1;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV1, // phiNew
			           this->tv(6), // src1 apply op on Im|alpha(C)>
			           direction);

			if (tmpV1.size() > 0)
				this->tvNonConst(6) = tmpV1;

			VectorWithOffsetType tmpV2;
			applyOneOp(loopNumber,
			           indexOfOperator,
			           site,
			           tmpV2,                            // phiNew
			           this->tv(7), // src1 apply op on Re|alpha(C)>
			           direction);

			if (tmpV2.size() > 0) {
				this->tvNonConst(7) = tmpV2;
				applied_ = true;
				PsimagLite::OstringStream msgg(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg = msgg();
				msg<<"PROD: Second Operator Applied";
				progress_.printline(msgg, std::cout);
			}
		}
	}

	void cocoon(SizeType site, ProgramGlobals::DirectionEnum direction) const
	{
		const AlgorithmEnumType algo = tstStruct_.algorithm();
		const bool isChevy = (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV);

		SizeType nineOrTenOrFifteen = (isChevy) ? 10 : 15;
		SizeType eightOrEleven = (isChevy) ? 8 : 11;

		if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			nineOrTenOrFifteen = 9;
			eightOrEleven = 8;
		}

		const ComplexOrRealType rr =
		        this->common().rixsCocoon(direction,site,nineOrTenOrFifteen,5,true);
		const ComplexOrRealType ri =
		        this->common().rixsCocoon(direction,site,nineOrTenOrFifteen,4,true);
		const ComplexOrRealType ir =
		        this->common().rixsCocoon(direction,site,eightOrEleven,5,true);
		const ComplexOrRealType ii =
		        this->common().rixsCocoon(direction,site,eightOrEleven,4,true);

		RealType time = 0.0;
		if (tstStruct_.algorithm() != TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			time = this->common().aoe().timeVectors().time();
		}
		std::cout<<site<<" "<<(ri-ir)<<" "<<time; // time here is the currentTime
		std::cout<<" <gs|A|P2> 1\n";   // 1 here is the "superdensity"
		std::cout<<site<<" "<<(rr+ii)<<" "<<time; // time here is the currentTime
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
			PsimagLite::OstringStream msgg2(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
			msg2<<"Norm of phi is zero\n";
			progress_.printline(msgg2, std::cout);
		}
	}

	void calcDynVectors(RealType Eg,
	                    ProgramGlobals::DirectionEnum direction,
	                    const VectorSizeType& block1)
	{
		if (!applied_ && appliedFirst_) {
			return;
		}

		if (!applied_) {
			return;
		}

		const AlgorithmEnumType algo = tstStruct_.algorithm();

		if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
			skeleton_.calcDynVectors(this->tv(6),
			                         this->tv(7),
			                         this->tvNonConst(8),
			                         this->tvNonConst(9));
			firstCall_ = false; // unused here but just in case
			return;
		}

		VectorSizeType indices;
		VectorSizeType indices2;
		//SizeType numberOfWeights = 0;

		if (algo == TargetParamsType::BaseType::AlgorithmEnum::CHEBYSHEV) {
			indices = {6, 8, 9};
			indices2 = {7, 10, 11};
			//numberOfWeights = 12;
		} else if (algo == TargetParamsType::BaseType::AlgorithmEnum::KRYLOVTIME){
			indices = {6, 8, 9, 10, 11};
			indices2 = {7, 12, 13, 14, 15};
			//numberOfWeights = 16;
		}

		//assert(numberOfWeights > 0);
		assert(indices.size() > 0 && indices2.size() > 0);
		calcVectors(indices, Eg, direction, block1, !firstCall_, false);
		calcVectors(indices2, Eg, direction, block1, !firstCall_, true);
		firstCall_ = false;
	}

	void calcVectors(const VectorSizeType& indices,
	                 RealType Eg,
	                 ProgramGlobals::DirectionEnum direction,
	                 const VectorSizeType& block1,
	                 bool wftOrAdvance,
	                 bool isLastCall)
	{
		bool allOperatorsApplied = (this->common().aoe().noStageIs(StageEnumType::DISABLED) &&
		                            this->common().aoe().noStageIs(StageEnumType::OPERATOR));

		const VectorWithOffsetType& v0 = this->tv(indices[0]);

		this->common().aoeNonConst().calcTimeVectors(indices,
		                                             Eg,
		                                             v0,
		                                             direction,
		                                             allOperatorsApplied,
		                                             wftOrAdvance, // wft and advance indices[0]
		                                             block1,
		                                             isLastCall);
	}

	void applyOneOp(SizeType loopNumber,
	                SizeType indexOfOperator,
	                SizeType site,
	                VectorWithOffsetType& dest,
	                const VectorWithOffsetType& src,
	                ProgramGlobals::DirectionEnum direction)
	{
		const AlgorithmEnumType algo = tstStruct_.algorithm();

		if (algo != TargetParamsType::BaseType::AlgorithmEnum::KRYLOV)
			this->common().aoeNonConst().applyOneOperator(loopNumber,
		                                                  indexOfOperator,
		                                                  site,
		                                                  dest, // phiNew
		                                                  src, // src1
		                                                  direction,
		                                                  *tstStruct2_);
		else
			this->common().aoeNonConst().applyOneOperator(loopNumber,
		                                                  indexOfOperator,
		                                                  site,
		                                                  dest, // phiNew
		                                                  src, // src1
		                                                  direction,
		                                                  tstStruct_);
	}

	void setWeights()
	{
		const SizeType n = this->numberOfTvs();
		weight_.resize(n);
		std::fill(weight_.begin(), weight_.end(), 0);
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			RealType norma = norm(this->tv(i));
			if (norma < 1e-6) continue;
			weight_[i] = 1/norma;
			sum += weight_[i];
		}

		gsWeightActual_ = 1 - sum;

		if (gsWeightActual_ >= gsWeight_) return; // <--- EARLY EXIT HERE

		assert(sum > 1e-6);
		RealType factor = (1 - gsWeight_)/sum;
		for (SizeType i = 0; i < n; ++i)
			weight_[i] *= factor;
		gsWeightActual_ = gsWeight_;
	}

	TargetParamsType tstStruct_;
	TargetParams2Type* tstStruct2_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	RealType gsWeightActual_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
	CorrectionVectorSkeletonType skeleton_;
	bool applied_;
	bool appliedFirst_;
	static bool firstCall_;
}; // class TargetingRixsDynamic

template<typename T1, typename T2>
bool TargetingRixsDynamic<T1, T2>::firstCall_ = true;

} // namespace
/*@}*/
#endif // TARGETING_RIXS_DYNAMIC_H
