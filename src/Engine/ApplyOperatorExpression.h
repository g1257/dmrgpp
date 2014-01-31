/*
Copyright (c) 2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
#ifndef APPLY_OP_EXPRESSION_H
#define APPLY_OP_EXPRESSION_H

#include "ProgressIndicator.h"
#include "ApplyOperatorLocal.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"

namespace Dmrg {

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
class ApplyOperatorExpression {

	typedef typename TargetHelperType::RealType RealType;
	typedef typename TargetHelperType::ModelType ModelType;
	typedef typename TargetHelperType::TargetParamsType TargetParamsType;
	typedef typename TargetHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetHelperType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef TimeVectorsBase<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsBaseType;
	typedef TimeVectorsKrylov<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsKrylovType;
	typedef TimeVectorsRungeKutta<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsRungeKuttaType;
	typedef TimeVectorsSuzukiTrotter<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsSuzukiTrotterType;

	static SizeType const PRODUCT = TargetParamsType::PRODUCT;
	static SizeType const SUM = TargetParamsType::SUM;

public:

	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType> ApplyOperatorType;
	typedef typename ApplyOperatorType::BorderEnum BorderEnumType;
	typedef typename TimeVectorsBaseType::PairType PairType;

	enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};

	ApplyOperatorExpression(const TargetHelperType& targetHelper,
	                        SizeType targets,
	                        SizeType indexNoAdvance)
	    : progress_("ApplyOperatorExpression"),
	      stage_(targetHelper.tstStruct().sites(),DISABLED),
	      E0_(0.0),
	      currentTime_(0.0),
	      indexNoAdvance_(indexNoAdvance),
	      applyOpLocal_(targetHelper.lrs()),
	      targetVectors_(targets),
	      timeVectorsBase_(0),
	      targetHelper_(targetHelper)
	{}

	~ApplyOperatorExpression()
	{
		if (timeVectorsBase_)
			delete timeVectorsBase_;
	}

	SizeType getPhi(VectorWithOffsetType& phiNew,
	                RealType Eg,
	                SizeType direction,
	                SizeType site,
	                SizeType loopNumber)
	{
		SizeType count =0;
		VectorWithOffsetType phiOld = psi_;
		VectorWithOffsetType vectorSum;

		SizeType max = targetHelper_.tstStruct().sites();
		if (noStageIs(DISABLED)) {
			max = 1;
			for (SizeType i=0;i<stage_.size();i++) {
				if (stage_[i]==OPERATOR) stage_[i] = WFT_NOADVANCE;
				if (stage_[i]==WFT_ADVANCE) stage_[i] = WFT_NOADVANCE;
			}
		}

		// Loop over each operator that needs to be applied
		// in turn to the g.s.
		for (SizeType i=0;i<max;i++) {

			SizeType count2 = evolve(i,phiNew,phiOld,Eg,direction,site,loopNumber,max-1);
			if (count2 == 0) continue;
			count += count2;

			if (targetHelper_.tstStruct().concatenation() == PRODUCT) {
				phiOld = phiNew;
			} else {
				if (stage_[i] == OPERATOR) vectorSum += phiNew;
				else vectorSum = phiNew;
			}
		}

		if (targetHelper_.tstStruct().concatenation() == SUM) phiNew = vectorSum;

		if (allStages(DISABLED)) E0_ = Eg;
		return count;
	}

	bool allStages(SizeType x) const
	{
		for (SizeType i=0;i<stage_.size();i++)
			if (stage_[i]!=x) return false;
		return true;
	}

	bool noStageIs(SizeType x) const
	{
		for (SizeType i=0;i<stage_.size();i++)
			if (stage_[i]==x) return false;
		return true;
	}

	const VectorSizeType& stage() const { return stage_; }

	void setAllStagesTo(SizeType x)
	{
		for (SizeType i=0;i<stage_.size();i++)
			stage_[i] = x;
	}

	const RealType& energy() const
	{
		return E0_;
	}

	const RealType& currentTime() const
	{
		return currentTime_;
	}

	const ApplyOperatorType& applyOpLocal() const
	{
		return applyOpLocal_;
	}

	const VectorSizeType& nonZeroQns() const
	{
		return nonZeroQns_;
	}

	VectorWithOffsetType& psi() // <--- FIXME
	{
		return psi_;
	}

	const VectorWithOffsetType& psi() const
	{
		return psi_;
	}

	const VectorVectorWithOffsetType& targetVectors() const
	{
		return targetVectors_;
	}

	VectorWithOffsetType& targetVectors(SizeType i) // <--- FIXME
	{
		return targetVectors_[i];
	}

	void targetVectorsResize(SizeType x)
	{
		targetVectors_.resize(x);
	}

	void initTimeVectors(const VectorRealType& times)
	{
		const LeftRightSuperType& lrs = targetHelper_.lrs();
		const ModelType& model = targetHelper_.model();
		const TargetParamsType& tstStruct = targetHelper_.tstStruct();
		const WaveFunctionTransfType& wft = targetHelper_.wft();

		PsimagLite::String s (__FILE__);
		s += " Unknown algorithm\n";

		switch (tstStruct.algorithm()) {
		case TargetParamsType::KRYLOV:
			timeVectorsBase_ = new TimeVectorsKrylovType(currentTime_,
			                                             tstStruct,
			                                             times,
			                                             targetVectors_,
			                                             model,
			                                             wft,
			                                             lrs,
			                                             E0_);
			break;
		case TargetParamsType::RUNGE_KUTTA:
			timeVectorsBase_ = new TimeVectorsRungeKuttaType(currentTime_,
			                                                 tstStruct,
			                                                 times,
			                                                 targetVectors_,
			                                                 model,
			                                                 wft,
			                                                 lrs,
			                                                 E0_);
			break;
		case TargetParamsType::SUZUKI_TROTTER:
			timeVectorsBase_ = new TimeVectorsSuzukiTrotterType(currentTime_,
			                                                    tstStruct,
			                                                    times,
			                                                    targetVectors_,
			                                                    model,
			                                                    wft,
			                                                    lrs,
			                                                    E0_,
			                                                    &nonZeroQns_);
			break;
		default:
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	void setTime(RealType t)
	{
		currentTime_ = t;
	}

	template<typename SomeSerializerType>
	void loadTargetVectors(SomeSerializerType& serializer)
	{
		for (SizeType i=0;i<targetVectors_.size();i++)
			targetVectors_[i] = serializer.vector(i);
	}

	void calcTimeVectors(const PairType& startEnd,
	                     RealType Eg,
	                     const VectorWithOffsetType& phi,
	                     SizeType direction,
	                     bool allOperatorsApplied,
	                     const PsimagLite::Vector<SizeType>::Type& block)
	{
		timeVectorsBase_->calcTimeVectors(startEnd,
		                                  Eg,
		                                  phi,
		                                  direction,
		                                  allOperatorsApplied,
		                                  block);
	}

private:

	void checkOrder(SizeType i) const
	{
		if (i==0) return;
		for (SizeType j=0;j<i;j++) {
			if (stage_[j] == DISABLED) {
				PsimagLite::String s ="TST:: Seeing dynamic site ";
				s += ttos(targetHelper_.tstStruct().sites(i));
				s =s + " before having seen";
				s = s + " site "+ttos(targetHelper_.tstStruct().sites(j));
				s = s +". Please order your dynamic sites in order of appearance.\n";
				throw PsimagLite::RuntimeError(s);
			}
		}
	}

	PsimagLite::String getStage(SizeType i) const
	{
		switch (stage_[i]) {
		case DISABLED:
			return "Disabled";
			break;
		case OPERATOR:
			return "Applying operator for the first time";
			break;
		case WFT_NOADVANCE:
			return "WFT_NOADVANCE";
			break;
		case WFT_ADVANCE:
			return "WFT_ADVANCE";
			break;
		}

		return "undefined";
	}

	void findElectronsOfOneSite(typename PsimagLite::Vector<SizeType>::Type& electrons,
	                            SizeType site) const
	{
		typename PsimagLite::Vector<SizeType>::Type block(1,site);
		typename ModelType::HilbertBasisType basis;
		typename PsimagLite::Vector<SizeType>::Type quantumNumbs;
		targetHelper_.model().setNaturalBasis(basis,quantumNumbs,block);
		targetHelper_.model().findElectrons(electrons,basis,site);
	}

	SizeType evolve(SizeType i,
	                VectorWithOffsetType& phiNew,
	                VectorWithOffsetType& phiOld,
	                RealType Eg,
	                SizeType direction,
	                SizeType site,
	                SizeType loopNumber,
	                SizeType lastI)
	{
		static SizeType timesWithoutAdvancement = 0;
		static bool firstSeeLeftCorner = false;
		SizeType advanceEach = targetHelper_.tstStruct().advanceEach();

		if (direction == ProgramGlobals::INFINITE) {
			E0_ = Eg;
			return 0;
		}

		if (targetHelper_.tstStruct().startingLoops().size()>0 &&
		    targetHelper_.tstStruct().startingLoops()[i]>loopNumber) return 0;

		if (site != targetHelper_.tstStruct().sites(i) && stage_[i]==DISABLED)
			return 0;

		if (site != targetHelper_.tstStruct().sites(i) && stage_[i]!=DISABLED && i>0)
			return 0;

		if (site == targetHelper_.tstStruct().sites(i) && stage_[i]==DISABLED) {
			stage_[i]=OPERATOR;
		} else {
			stage_[i]=WFT_NOADVANCE;
		}

		if (stage_[i] == OPERATOR) checkOrder(i);

		if (advanceEach > 0 && timesWithoutAdvancement >= advanceEach) {
			stage_[i] = WFT_ADVANCE;
			if (i==lastI) {
				currentTime_ += targetHelper_.tstStruct().tau();
				timesWithoutAdvancement=1;
			}
		} else {
			if (i==lastI && stage_[i]==WFT_NOADVANCE && firstSeeLeftCorner)
				timesWithoutAdvancement++;
		}

		if (!firstSeeLeftCorner && i==lastI && stage_[i]==WFT_NOADVANCE && site==1)
			firstSeeLeftCorner=true;

		PsimagLite::OstringStream msg2;
		msg2<<"Steps without advance: "<<timesWithoutAdvancement;
		msg2<<" site="<<site<<" currenTime="<<currentTime_;
		if (timesWithoutAdvancement>0) progress_.printline(msg2,std::cout);

		PsimagLite::OstringStream msg;
		msg<<"Evolving, stage="<<getStage(i);
		msg<<" site="<<site<<" loopNumber="<<loopNumber;
		msg<<" Eg="<<Eg;
		progress_.printline(msg,std::cout);

		// phi = A|psi>
		computePhi(i,site,phiNew,phiOld,direction);

		return 1;
	}

	void computePhi(SizeType i,
	                SizeType site,
	                VectorWithOffsetType& phiNew,
	                VectorWithOffsetType& phiOld,
	                SizeType systemOrEnviron)
	{
		SizeType numberOfSites = targetHelper_.lrs().super().block().size();
		SizeType advanceEach = targetHelper_.tstStruct().advanceEach();

		if (stage_[i]==OPERATOR) {

			BorderEnumType corner = (targetHelper_.tstStruct().sites(i)==0 ||
			                         targetHelper_.tstStruct().sites(i)==numberOfSites -1) ?
			            ApplyOperatorType::BORDER_YES : ApplyOperatorType::BORDER_NO;

			PsimagLite::OstringStream msg;
			msg<<"I'm applying a local operator now";
			progress_.printline(msg,std::cout);
			typename PsimagLite::Vector<SizeType>::Type electrons;
			findElectronsOfOneSite(electrons,site);
			FermionSign fs(targetHelper_.lrs().left(),electrons);
			applyOpLocal_(phiNew,phiOld,targetHelper_.tstStruct().aOperators()[i],
			              fs,systemOrEnviron,corner);
			RealType norma = std::norm(phiNew);

			if (norma<1e-6) {
				PsimagLite::OstringStream msg2;
				msg2<<"Norm of phi is zero\n";
				progress_.printline(msg2,std::cout);
			}

			if (targetHelper_.tstStruct().useQns()) setQuantumNumbers(phiNew);

		} else if (stage_[i] >= WFT_NOADVANCE) {

			SizeType advance = indexNoAdvance_;
			if (advanceEach > 0 && stage_[i] == WFT_ADVANCE) {
				SizeType timeSteps = targetHelper_.tstStruct().timeSteps();
				advance = (timeSteps > 0) ? timeSteps - 1 : 0;
				timeVectorsBase_->timeHasAdvanced();
			}

			if (site==0 || site==numberOfSites -1)  {
				// don't wft since we did it before
				phiNew = targetVectors_[advance];
				return;
			}

			PsimagLite::OstringStream msg;
			msg<<"I'm calling the WFT now";
			progress_.printline(msg,std::cout);

			if (targetHelper_.tstStruct().aOperators().size() == 1) {
				guessPhiSectors(phiNew,i,systemOrEnviron,site);
			} else {
				if (targetHelper_.tstStruct().useQns())
					phiNew.populateFromQns(nonZeroQns_,targetHelper_.lrs().super());
				else
					phiNew.populateSectors(targetHelper_.lrs().super());
			}

			// OK, now that we got the partition number right, let's wft:
			VectorSizeType nk(1,targetHelper_.model().hilbertSize(site));
			targetHelper_.wft().setInitialVector(phiNew,
			                                     targetVectors_[advance],
			                                     targetHelper_.lrs(),
			                                     nk);
			phiNew.collapseSectors();

		} else {
			throw PsimagLite::RuntimeError("computePhi\n");
		}
	}

	void guessPhiSectors(VectorWithOffsetType& phi,
	                     SizeType i,
	                     SizeType systemOrEnviron,
	                     SizeType site)
	{
		VectorSizeType electrons;
		targetHelper_.model().findElectronsOfOneSite(electrons,site);
		FermionSign fs(targetHelper_.lrs().left(),electrons);
		if (allStages(WFT_NOADVANCE)) {
			VectorWithOffsetType tmpVector = psi_;
			for (SizeType j=0;j<targetHelper_.tstStruct().aOperators().size();j++) {
				applyOpLocal_(phi,
				              tmpVector,
				              targetHelper_.tstStruct().aOperators()[j],
				              fs,
				              systemOrEnviron,
				              ApplyOperatorType::BORDER_NO);
				tmpVector = phi;
			}
			return;
		}
		applyOpLocal_(phi,
		              psi_,
		              targetHelper_.tstStruct().aOperators()[i],
		              fs,
		              systemOrEnviron,
		              ApplyOperatorType::BORDER_NO);
	}

	void setQuantumNumbers(const VectorWithOffsetType& v)
	{
		nonZeroQns_.clear();
		assert(v.size()==targetHelper_.lrs().super().size());
		for (SizeType i=0;i<v.sectors();i++) {
			SizeType i0 = v.sector(i);
			SizeType state = v.offset(i0);
			SizeType qn = targetHelper_.lrs().super().qn(state);
			nonZeroQns_.push_back(qn);
		}
	}

	PsimagLite::ProgressIndicator progress_;
	VectorSizeType stage_;
	RealType E0_;
	RealType currentTime_;
	SizeType indexNoAdvance_;
	ApplyOperatorType applyOpLocal_;
	VectorSizeType nonZeroQns_;
	VectorWithOffsetType psi_;
	typename PsimagLite::Vector<VectorWithOffsetType>::Type targetVectors_;
	TimeVectorsBaseType* timeVectorsBase_;
	const TargetHelperType& targetHelper_;
};

} // namespace Dmrg

#endif // APPLY_OP_EXPRESSION_H

