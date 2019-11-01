/*
Copyright (c) 2014-2017-2018-2019, UT-Battelle, LLC
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
#ifndef APPLY_OP_EXPRESSION_H
#define APPLY_OP_EXPRESSION_H

#include "ProgressIndicator.h"
#include "ApplyOperatorLocal.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsChebyshev.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include "Io/IoSelector.h"
#include "TargetParamsBase.h"
#include "StageEnum.h"
#include "MultiSiteExpressionHelper.h"
#include "CorrelationsSkeleton.h"

namespace Dmrg {

template<typename TargetHelperType,
         typename VectorWithOffsetType,
         typename LanczosSolverType>
class ApplyOperatorExpression {

public:

	typedef typename TargetHelperType::RealType RealType;
	typedef typename TargetHelperType::ModelType ModelType;
	typedef typename ModelType::ParametersType ParametersType;
	typedef typename ParametersType::OptionsType OptionsType;
	typedef TargetParamsBase<ModelType> TargetParamsType;
	typedef typename TargetHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename TargetHelperType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef TimeVectorsBase<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsBaseType;
	typedef TimeVectorsKrylov<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsKrylovType;
	typedef TimeVectorsChebyshev<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsChebyshevType;
	typedef TimeVectorsRungeKutta<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsRungeKuttaType;
	typedef TimeVectorsSuzukiTrotter<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsSuzukiTrotterType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType> ApplyOperatorType;
	typedef typename ApplyOperatorType::BorderEnum BorderEnumType;
	typedef typename TimeVectorsBaseType::PairType PairType;
	typedef StageEnum StageEnumType;
	typedef typename PsimagLite::Vector<StageEnum>::Type VectorStageEnumType;
	typedef MultiSiteExpressionHelper<LeftRightSuperType, VectorWithOffsetType>
	MultiSiteExpressionHelperType;
	typedef typename MultiSiteExpressionHelperType::DmrgSerializerType DmrgSerializerType;
	typedef CorrelationsSkeleton<MultiSiteExpressionHelperType, ModelType> CorrelationsSkeletonType;
	typedef typename TimeVectorsBaseType::WftHelperType WftHelperType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename PsimagLite::Vector<typename
	PsimagLite::Vector<VectorWithOffsetType*>::Type>::Type VectorVectorVectorWithOffsetType;

	ApplyOperatorExpression(const TargetHelperType& targetHelper,
	                        SizeType indexNoAdvance)
	    : progress_("ApplyOperatorExpression"),
	      targetHelper_(targetHelper),
	      E0_(0.0),
	      currentTimeStep_(0),
	      indexNoAdvance_(indexNoAdvance),
	      applyOpLocal_(targetHelper.lrs(), targetHelper.withLegacyBugs()),
	      targetVectors_(0),
	      timeVectorsBase_(0),
	      wftHelper_(targetHelper.model(), targetHelper.lrs(), targetHelper.wft()),
	      multiSiteExprHelper_(targetHelper_.model().geometry().numberOfSites() - 2),
	      correlationsSkel_(multiSiteExprHelper_, false)
	{}

	~ApplyOperatorExpression()
	{
		delete timeVectorsBase_;
		timeVectorsBase_ = 0;
		clearPsi();
	}

	virtual void initPsi(SizeType nexcited, SizeType nsectors)
	{
		if (psi_.size() > 0) {
			if (psi_.size() != nexcited)
				err("AOE::initPsi(): nexcited cannot change during run\n");
			bool flag = true;
			for (SizeType i = 0; i < nexcited; ++i)
				if (psi_[i].size() != nsectors)
					flag = false;

			if (flag) return;

			std::cerr<<"AOE::initPsi(): WARNING sectors changed during run\n";
			std::cout<<"AOE::initPsi(): WARNING sectors changed during run\n";
			for (SizeType i = 0; i < nexcited; ++i) {
				for (SizeType j = 0; j < psi_[i].size(); ++j) {
					delete psi_[i][j];
					psi_[i][j] = nullptr;
				}
			}
		}

		psi_.resize(nexcited);

		for (SizeType i = 0; i < nexcited; ++i)
			psi_[i].resize(nsectors, nullptr);
	}


	void postCtor(SizeType tstSites)
	{
		if (stage_.size() != 0)
			throw PsimagLite::RuntimeError("ApplyOperatorExpression: Internal Error\n");

		stage_.resize((tstSites == 0) ? 1 : tstSites, StageEnum::DISABLED);
	}

	SizeType getPhi(VectorWithOffsetType* phiNew,
	                RealType Eg,
	                ProgramGlobals::DirectionEnum direction,
	                SizeType site,
	                SizeType loopNumber,
	                const TargetParamsType& tstStruct)
	{
		SizeType count =0;
		VectorWithOffsetType phiOld = ensureOnlyOnePsi("getPhi");
		VectorWithOffsetType vectorSum;

		SizeType max = tstStruct.sites();
		if (noStageIs(StageEnum::DISABLED)) {
			max = 1;
			for (SizeType i=0;i<stage_.size();i++) {
				if (stage_[i] == StageEnum::OPERATOR)
					stage_[i] = StageEnum::WFT_NOADVANCE;
				if (stage_[i] == StageEnum::WFT_ADVANCE)
					stage_[i] = StageEnum::WFT_NOADVANCE;
			}
		}

		// Loop over each operator that needs to be applied
		// in turn to the g.s.
		for (SizeType i = 0; i < max; ++i) {

			SizeType count2 = evolve(i, Eg, direction, site, loopNumber, max - 1, tstStruct);

			if (count2 == 0) continue;

			// phi = A|psi>
			if (phiNew)
				computePhi(i, site, *phiNew, phiOld, direction, tstStruct);

			count += count2;

			if (!phiNew) continue;

			if (tstStruct.concatenation() == TargetParamsType::ConcatEnum::PRODUCT) {
				phiOld = *phiNew;
			} else {
				if (stage_[i] == StageEnum::OPERATOR) vectorSum += *phiNew;
				else vectorSum = *phiNew;
			}
		}

		if (phiNew && tstStruct.concatenation() == TargetParamsType::ConcatEnum::SUM)
			*phiNew = vectorSum;

		if (allStages(StageEnum::DISABLED)) E0_ = Eg;

		if (noStageIs(StageEnum::DISABLED)) {
			PsimagLite::OstringStream msg;
			msg<<"EnergyForExp was ";
			if (tstStruct.isEnergyForExp()) {
				E0_ = tstStruct.energyForExp();
				msg<<"provided explicitly, ";
			} else {
				msg<<"not provided explicitly, ";
			}

			msg<<" value= "<<E0_;
			progress_.printline(msg, std::cout);
		}

		return count;
	}

	bool allStages(StageEnum x) const
	{
		for (SizeType i=0;i<stage_.size();i++)
			if (stage_[i] != x) return false;
		return true;
	}

	bool noStageIs(StageEnum x) const
	{
		for (SizeType i=0;i<stage_.size();i++)
			if (stage_[i] == x) return false;
		return true;
	}

	const VectorStageEnumType& stages() const { return stage_; }

	void setStage(SizeType ind, StageEnum x)
	{
		assert(ind < stage_.size());
		stage_[ind] = x;
	}

	const RealType& energy() const
	{
		return E0_;
	}

	const ApplyOperatorType& applyOpLocal() const
	{
		return applyOpLocal_;
	}

	const VectorVectorVectorWithOffsetType& psi() const
	{
		return psi_;
	}

	void setOnlyOnePsi(const VectorWithOffsetType& v)
	{
		std::cerr<<"Setting only one psi\n";
		std::cout<<"Setting only one psi\n";

		clearPsi();

		psi_.resize(1);
		psi_[0].resize(1);
		psi_[0][0] = new VectorWithOffsetType(v);
	}

	template<typename SomeBasisType>
	void setPsi(SizeType excitedIndex,
	            SizeType sectorIndex,
	            VectorType& v,
	            const SomeBasisType& basis,
	            const VectorSizeType& sectors)
	{
		if (psi_.size() == 0) {
			psi_.resize(targetHelper_.model().params().excited.size());
			for (SizeType i = 0; i < psi_.size(); ++i)
				psi_[i].resize(sectors.size());
		}

		assert(sectorIndex < sectors.size());
		const SizeType sector = sectors[sectorIndex];

		assert(excitedIndex < psi_.size());
		assert(sectorIndex < psi_[excitedIndex].size());

		if (!psi_[excitedIndex][sectorIndex])
			psi_[excitedIndex][sectorIndex] = new VectorWithOffsetType;
		psi_[excitedIndex][sectorIndex]->set(v, sector, basis);
	}

	void writePsi(PsimagLite::IoSelector::Out& io, PsimagLite::String prefix) const
	{
		const SizeType nexcited = psi_.size();
		io.createGroup(prefix + "/PSI");
		io.write(nexcited, prefix + "/PSI/Size");
		for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {
			const SizeType nsectors = psi_[excitedIndex].size();
			PsimagLite::String label = prefix + "/PSI/" + ttos(excitedIndex);
			io.createGroup(label);
			io.write(nsectors, label + "/Size");
			for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex)
				psi_[excitedIndex][sectorIndex]->write(io, label + "/" + ttos(sectorIndex));
		}
	}

	void readPsi(PsimagLite::IoSelector::In& io, PsimagLite::String prefix)
	{
		SizeType nexcited = 0;
		io.read(nexcited, prefix + "/PSI/Size");
		clearPsi();
		psi_.resize(nexcited);
		for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {

			PsimagLite::String label = prefix + "/PSI/" + ttos(excitedIndex) + "/";
			SizeType nsectors = 0;
			io.read(nsectors, label + "Size");
			psi_[excitedIndex].resize(nsectors);
			for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
				psi_[excitedIndex][sectorIndex] = new VectorWithOffsetType;
				psi_[excitedIndex][sectorIndex]->read(io, label + ttos(sectorIndex));
			}
		}
	}

	const VectorVectorWithOffsetType& targetVectors() const
	{
		return targetVectors_;
	}

	VectorWithOffsetType& targetVectors(SizeType i) // <--- FIXME
	{
		assert(i < targetVectors_.size());
		return targetVectors_[i];
	}

	const VectorWithOffsetType& targetVectors(SizeType i) const
	{
		assert(i < targetVectors_.size());
		return targetVectors_[i];
	}

	void targetVectorsResize(SizeType x)
	{
		if (x == 0) targetVectors_.clear();
		else targetVectors_.resize(x);
	}

	void initTimeVectors(const TargetParamsType& tstStruct,
	                     const VectorRealType& times,
	                     InputValidatorType& ioIn)
	{
		const LeftRightSuperType& lrs = targetHelper_.lrs();
		const ModelType& model = targetHelper_.model();
		const WaveFunctionTransfType& wft = targetHelper_.wft();

		PsimagLite::String s (__FILE__);
		s += " Unknown algorithm\n";

		switch (tstStruct.algorithm()) {
		case TargetParamsType::AlgorithmEnum::KRYLOV:
			timeVectorsBase_ = new TimeVectorsKrylovType(currentTimeStep_,
			                                             tstStruct,
			                                             times,
			                                             targetVectors_,
			                                             model,
			                                             wft,
			                                             lrs,
			                                             ioIn);
			break;
		case TargetParamsType::AlgorithmEnum::CHEBYSHEV:
			timeVectorsBase_ = new TimeVectorsChebyshevType(currentTimeStep_,
			                                                tstStruct,
			                                                times,
			                                                targetVectors_,
			                                                model,
			                                                wft,
			                                                lrs);
			break;
		case TargetParamsType::AlgorithmEnum::RUNGE_KUTTA:
			timeVectorsBase_ = new TimeVectorsRungeKuttaType(currentTimeStep_,
			                                                 tstStruct,
			                                                 times,
			                                                 targetVectors_,
			                                                 model,
			                                                 wft,
			                                                 lrs);
			break;
		case TargetParamsType::AlgorithmEnum::SUZUKI_TROTTER:
			timeVectorsBase_ = new TimeVectorsSuzukiTrotterType(currentTimeStep_,
			                                                    tstStruct,
			                                                    times,
			                                                    targetVectors_,
			                                                    model,
			                                                    wft,
			                                                    lrs);
			break;
		default:
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	RealType time() const
	{
		return (timeVectorsBase_) ? timeVectorsBase_->time() : 0;
	}

	SizeType currentTimeStep() const { return currentTimeStep_; }

	void setCurrentTimeStep(SizeType t)
	{
		currentTimeStep_ = t;
	}

	void timeHasAdvanced() { timeVectorsBase_->timeHasAdvanced(); }

	void loadEnergy(PsimagLite::IoSelector::In& io,
	                PsimagLite::String label)
	{
		io.readLastVectorEntry(E0_, label);
	}

	void multiplyTimeVector(SizeType i,RealType factor)
	{
		targetVectors_[i] = factor*targetVectors_[i];
	}

	void calcTimeVectors(const PsimagLite::Vector<SizeType>::Type& indices,
	                     RealType Eg,
	                     const VectorWithOffsetType& phi,
	                     ProgramGlobals::DirectionEnum direction,
	                     bool allOperatorsApplied,
	                     bool wftAndAdvanceIfNeeded,
	                     const PsimagLite::Vector<SizeType>::Type& block,
	                     bool isLastCall)
	{
		typename TimeVectorsBaseType::ExtraData extra(direction,
		                                              allOperatorsApplied,
		                                              wftAndAdvanceIfNeeded,
		                                              block,
		                                              isLastCall);
		timeVectorsBase_->calcTimeVectors(indices,
		                                  Eg,
		                                  phi,
		                                  extra);
	}

	void applyOneOperator(SizeType loopNumber,
	                      SizeType indexOfOperator,
	                      SizeType site,
	                      VectorWithOffsetType& phiNew,
	                      const VectorWithOffsetType& psiSrc,
	                      const ProgramGlobals::DirectionEnum systemOrEnviron,
	                      const TargetParamsType& tstStruct)
	{
		if (tstStruct.startingLoops().size()>0 &&
		        tstStruct.startingLoops()[indexOfOperator]>loopNumber)
			return;

		const bool hasBeenApplied = (phiNew.size() > 0);
		if (hasBeenApplied) return;

		VectorWithOffsetType phiOld = psiSrc;
		SizeType numberOfSites = targetHelper_.lrs().super().block().size();

		BorderEnumType corner = (site == 0 || site == numberOfSites -1) ?
		            ApplyOperatorType::BORDER_YES : ApplyOperatorType::BORDER_NO;

		PsimagLite::OstringStream msg;
		msg<<"I'm applying a local operator now";
		progress_.printline(msg,std::cout);
		typename PsimagLite::Vector<bool>::Type signs;
		targetHelper_.model().findOddElectronsOfOneSite(signs, site);
		FermionSign fs(targetHelper_.lrs().left(), signs);
		applyOpLocal_(phiNew,
		              phiOld,
		              tstStruct.aOperators()[indexOfOperator],
		              fs,
		              systemOrEnviron,
		              corner);

		RealType norma = norm(phiNew);
		if (norma<1e-6) {
			PsimagLite::OstringStream msg2;
			msg2<<"Norm of phi is zero\n";
			progress_.printline(msg2,std::cout);
		}
	}

	void multiSitePush(DmrgSerializerType const* ds) const
	{
		const VectorWithOffsetType& psi00 = ensureOnlyOnePsi("multiSitePush");
		multiSiteExprHelper_.push(ds, psi00);
	}

	void wftSome(SizeType site, SizeType begin, SizeType end)
	{
		VectorVectorWithOffsetType& tvs = const_cast<VectorVectorWithOffsetType&>(targetVectors_);
		wftHelper_.wftSome(tvs, site, begin, end);
	}

	const VectorWithOffsetType& ensureOnlyOnePsi(PsimagLite::String func) const
	{
		if (psi_.size() != 1)
			err("ensureOnlyOnePsi failed (more than one excited); called from" +
			    func + "\n");

		if (psi_[0].size() != 1)
			err("ensureOnlyOnePsi failed (more than one sector); called from" +
			    func + "\n");

		return *(psi_[0][0]);
	}

private:

	void checkOrder(SizeType i, const TargetParamsType& tstStruct) const
	{
		if (i==0) return;
		for (SizeType j=0;j<i;j++) {
			if (stage_[j] == StageEnum::DISABLED) {
				PsimagLite::String s ="TST:: Seeing dynamic site ";
				s += ttos(tstStruct.sites(i));
				s =s + " before having seen";
				s = s + " site "+ttos(tstStruct.sites(j));
				s = s +". Please order your dynamic sites in order of appearance.\n";
				throw PsimagLite::RuntimeError(s);
			}
		}
	}

	PsimagLite::String stageToString(SizeType i) const
	{
		switch (stage_[i]) {
		case StageEnum::DISABLED:
			return "Disabled";
		case StageEnum::OPERATOR:
			return "Applying operator for the first time";
		case StageEnum::WFT_NOADVANCE:
			return "WFT_NOADVANCE";
		case StageEnum::WFT_ADVANCE:
			return "WFT_ADVANCE";
		case StageEnum::COLLAPSE:
			return "COLLAPSE";
		}

		return "undefined";
	}

	SizeType evolve(SizeType i,
	                RealType Eg,
	                ProgramGlobals::DirectionEnum direction,
	                SizeType site,
	                SizeType loopNumber,
	                SizeType lastI,
	                const TargetParamsType& tstStruct)
	{
		static SizeType timesWithoutAdvancement = 0;
		static bool firstSeeLeftCorner = false;
		SizeType advanceEach = tstStruct.advanceEach();

		if (direction == ProgramGlobals::DirectionEnum::INFINITE) {
			E0_ = Eg;
			return 0;
		}

		if (tstStruct.startingLoops().size()>0 &&
		        tstStruct.startingLoops()[i]>loopNumber) return 0;

		if (site != tstStruct.sites(i) && stage_[i] == StageEnum::DISABLED)
			return 0;

		if (site != tstStruct.sites(i) && stage_[i] != StageEnum::DISABLED && i>0)
			return 0;

		if (site == tstStruct.sites(i) && stage_[i] == StageEnum::DISABLED) {
			stage_[i] = StageEnum::OPERATOR;
		} else {
			stage_[i] = StageEnum::WFT_NOADVANCE;
		}

		if (stage_[i] == StageEnum::OPERATOR) checkOrder(i, tstStruct);

		const OptionsType& options = targetHelper_.model().params().options;
		bool advanceOnlyAtBorder = (options.find("advanceUnrestricted") ==
		                            PsimagLite::String::npos);
		SizeType sites = targetHelper_.model().geometry().numberOfSites();
		bool weAreAtBorder = (site < 2 || site >= sites-2);
		bool dontAdvance = (advanceOnlyAtBorder & !weAreAtBorder);

		if (advanceEach > 0 && timesWithoutAdvancement >= advanceEach && !dontAdvance) {
			stage_[i] = StageEnum::WFT_ADVANCE;
			if (i == lastI) {
				++currentTimeStep_;
				timesWithoutAdvancement=1;
				timeVectorsBase_->timeHasAdvanced();
			}
		} else {
			if (i == lastI &&
			        stage_[i] == StageEnum::WFT_NOADVANCE &&
			        firstSeeLeftCorner)
				timesWithoutAdvancement++;
		}

		if (!firstSeeLeftCorner &&
		        i==lastI &&
		        stage_[i] == StageEnum::WFT_NOADVANCE &&
		        site==1)
			firstSeeLeftCorner=true;

		PsimagLite::OstringStream msg2;
		msg2<<"Steps without advance: "<<timesWithoutAdvancement;
		msg2<<" site="<<site<<" currenTime="<<time();
		if (timesWithoutAdvancement>0) progress_.printline(msg2,std::cout);

		PsimagLite::OstringStream msg;
		msg<<"Evolving, stage="<<stageToString(i);
		msg<<" site="<<site<<" loopNumber="<<loopNumber;
		msg<<" Eg="<<Eg;
		progress_.printline(msg,std::cout);

		return 1;
	}

	void computePhi(SizeType i,
	                SizeType site,
	                VectorWithOffsetType& phiNew,
	                VectorWithOffsetType& phiOld,
	                const ProgramGlobals::DirectionEnum systemOrEnviron,
	                const TargetParamsType& tstStruct) const
	{
		SizeType numberOfSites = targetHelper_.lrs().super().block().size();

		if (stage_[i] == StageEnum::OPERATOR) {

			BorderEnumType corner = (tstStruct.sites(i)==0 ||
			                         tstStruct.sites(i)==numberOfSites -1) ?
			            ApplyOperatorType::BORDER_YES : ApplyOperatorType::BORDER_NO;

			PsimagLite::OstringStream msg;
			msg<<"I'm applying a local operator now";
			progress_.printline(msg,std::cout);
			typename PsimagLite::Vector<bool>::Type signs;
			targetHelper_.model().findOddElectronsOfOneSite(signs,site);
			FermionSign fs(targetHelper_.lrs().left(), signs);
			applyOpLocal_(phiNew,
			              phiOld,
			              tstStruct.aOperators()[i],
			              fs,
			              systemOrEnviron,corner);
			RealType norma = norm(phiNew);

			if (norma<1e-6) {
				PsimagLite::OstringStream msg2;
				msg2<<"Norm of phi is zero\n";
				progress_.printline(msg2,std::cout);
			}
		} else if (stage_[i] == StageEnum::WFT_NOADVANCE ||
		           stage_[i] == StageEnum::WFT_ADVANCE) {
			const SizeType advanceEach = tstStruct.advanceEach();
			SizeType advance = indexNoAdvance_;

			if (advanceEach > 0 && stage_[i] == StageEnum::WFT_ADVANCE) {
				SizeType timeSteps = tstStruct.timeSteps();
				advance = (timeSteps > 0) ? timeSteps - 1 : 0;
			}

			if (targetVectors_.size() <= advance) {
				PsimagLite::String s(__FILE__);
				s += ": TargetVectors.size()" + ttos(targetVectors_.size());
				s += " but advance=" + ttos(advance) + "\n";
				throw PsimagLite::RuntimeError(s);
			}

			const VectorWithOffsetType& src = targetVectors_[advance];

			if (src.size() == 0) {
				PsimagLite::String s(__FILE__);
				s += ": TargetVectors[" + ttos(advance) + "].size()==0\n";
				throw PsimagLite::RuntimeError(s);
			}

			if (site==0 || site==numberOfSites -1)  {
				// don't wft since we did it before
				assert(advance < targetVectors_.size());
				phiNew = src;
				return;
			}

			PsimagLite::OstringStream msg;
			msg<<"I'm calling the WFT now";

			wftHelper_.wftOneVector(phiNew, src, site);
		} else {
			throw PsimagLite::RuntimeError("computePhi\n");
		}
	}

	void clearPsi()
	{
		const SizeType nexcited = psi_.size();
		for (SizeType i =0; i < nexcited; ++i) {
			const SizeType sectors = psi_[i].size();
			for (SizeType j = 0; j < sectors; ++j) {
				delete psi_[i][j];
				psi_[i][j] = 0;
			}
		}
	}

	ApplyOperatorExpression(const ApplyOperatorExpression&) = delete;

	ApplyOperatorExpression& operator=(const ApplyOperatorExpression&) = delete;

	PsimagLite::ProgressIndicator progress_;
	const TargetHelperType& targetHelper_;
	VectorStageEnumType stage_;
	RealType E0_;
	SizeType currentTimeStep_;
	SizeType indexNoAdvance_;
	ApplyOperatorType applyOpLocal_;
	VectorVectorVectorWithOffsetType psi_;
	VectorVectorWithOffsetType targetVectors_;
	TimeVectorsBaseType* timeVectorsBase_;
	WftHelperType wftHelper_;
	mutable MultiSiteExpressionHelperType multiSiteExprHelper_;
	CorrelationsSkeletonType correlationsSkel_;
};

} // namespace Dmrg

#endif // APPLY_OP_EXPRESSION_H

