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
	typedef typename PsimagLite::Vector<VectorWithOffsetType*>::Type VectorVectorWithOffsetType;
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
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Vector<VectorVectorType>::Type VectorVectorVectorType;

	ApplyOperatorExpression(const TargetHelperType& targetHelper,
	                        SizeType indexNoAdvance)
	    : progress_("ApplyOperatorExpression"),
	      targetHelper_(targetHelper),
	      E0_(0.0),
	      indexNoAdvance_(indexNoAdvance),
	      applyOpLocal_(targetHelper.lrs(), targetHelper.withLegacyBugs()),
	      targetVectors_(0),
	      timeVectorsBase_(new TimeVectorsBaseType(targetHelper_.model(),
	                                               targetHelper_.lrs(),
	                                               targetHelper_.wft(),
	                                               "base")),
	      wftHelper_(targetHelper.model(), targetHelper.lrs(), targetHelper.wft()),
	      multiSiteExprHelper_(targetHelper_.model().superGeometry().numberOfSites() - 2),
	      correlationsSkel_(multiSiteExprHelper_, targetHelper.model(), false)
	{
		timesWithoutAdvancement_ = 0;
		firstSeeLeftCorner_ = false;
	}

	~ApplyOperatorExpression()
	{
		delete timeVectorsBase_;
		timeVectorsBase_ = 0;
		clearPsi();

		const SizeType n = targetVectors_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete targetVectors_[i];
			targetVectors_[i] = nullptr;
		}
	}

	void initPsi(SizeType nsectors, SizeType nexcited)
	{
		if (psi_.size() > 0) {
			bool flag = true;

			if (psi_.size() != nsectors) {
				flag = false;
			} else {
				for (SizeType i = 0; i < nsectors; ++i)
					if (psi_[i].size() != nexcited)
						flag = false;
			}

			if (flag) return;

			std::cerr<<"AOE::initPsi(): WARNING sectors/nexcited changed during run\n";
			std::cout<<"AOE::initPsi(): WARNING sectors/nexcited changed during run\n";
			clearPsi();
		}

		psi_.resize(nsectors);

		for (SizeType i = 0; i < nsectors; ++i)
			psi_[i].resize(nexcited, nullptr);
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
		SizeType count = 0;
		const SizeType nsectors = psi_.size();
		const SizeType sectorIndex = tstStruct.sectorIndex();
		if (sectorIndex >= nsectors)
			err("getPhi: sectors=" + ttos(nsectors) + " <= " + ttos(sectorIndex) + "\n");

		const SizeType nlevels = psi_[sectorIndex].size();
		const SizeType levelIndex = tstStruct.levelIndex();
		if (levelIndex >= nlevels)
			err("getPhi: levels=" + ttos(nlevels) + " <= " + ttos(levelIndex) + "\n");

		// deep copy needed
		VectorWithOffsetType phiOld = *(psi_[sectorIndex][levelIndex]);
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
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"EnergyForExp was ";
			if (tstStruct.isEnergyForExp()) {
				E0_ = tstStruct.energyForExp();
				msg<<"provided explicitly, ";
			} else {
				msg<<"not provided explicitly, ";
			}

			msg<<" value= "<<E0_;
			progress_.printline(msgg, std::cout);
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

	const RealType& energy() const { return E0_; }

	const ApplyOperatorType& applyOpLocal() const { return applyOpLocal_; }

	const VectorVectorVectorWithOffsetType& psiConst() const { return psi_; }

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
	void setPsi(VectorVectorVectorType& inV,
	            const VectorSizeType& sectors,
	            const SomeBasisType& someBasis,
	            SizeType nexcited)
	{
		const SizeType nsectors = sectors.size();

		assert(nsectors > 0);

		assert(nexcited > 0);

		if (nsectors != inV.size())
			err("FATAL: inV.size == " + ttos(inV.size()) + " but params.excited.size= "
			    + ttos(nsectors) + "\n");

		this->initPsi(nsectors, nexcited);

		for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
			if (inV[sectorIndex].size() != nexcited)
				err("Expected inV[" + ttos(sectorIndex) + "].size == " +
				    ttos(nexcited) + " but found " + ttos(inV[sectorIndex].size()) +
				    " instead\n");

			for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {

				if (!psi_[sectorIndex][excitedIndex]) {
					psi_[sectorIndex][excitedIndex] = new VectorWithOffsetType;
				}

				psi_[sectorIndex][excitedIndex]->set(inV[sectorIndex][excitedIndex],
				                                     sectorIndex,
				                                     someBasis);
			}
		}
	}

	void writePsi(PsimagLite::IoSelector::Out& io, PsimagLite::String prefix) const
	{
		const SizeType nsectors = psi_.size();
		io.createGroup(prefix + "/PSI");
		io.write(nsectors, prefix + "/PSI/Size");
		for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
			const SizeType nexcited = psi_[sectorIndex].size();
			PsimagLite::String label = prefix + "/PSI/" + ttos(sectorIndex);
			io.createGroup(label);
			io.write(nexcited, label + "/Size");
			for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex)
				psi_[sectorIndex][excitedIndex]->write(io, label + "/" + ttos(excitedIndex));
		}
	}

	void readPsi(PsimagLite::IoSelector::In& io, PsimagLite::String prefix)
	{
		clearPsi();
		SizeType nsectors = 0;
		try {
			io.read(nsectors, prefix + "/PSI/Size");
		} catch (...) {
			psi_.resize(1);
			psi_[0].resize(1);
			psi_[0][0] = new VectorWithOffsetType;
			psi_[0][0]->read(io, prefix + "PSI");
			return;
		}

		psi_.resize(nsectors);

		for (SizeType sectorIndex = 0; sectorIndex < nsectors; ++sectorIndex) {
			PsimagLite::String label = prefix + "/PSI/" + ttos(sectorIndex) + "/";
			SizeType nexcited = 0;
			io.read(nexcited, label + "Size");
			psi_[sectorIndex].resize(nexcited);
			for (SizeType excitedIndex = 0; excitedIndex < nexcited; ++excitedIndex) {
				psi_[sectorIndex][excitedIndex] = new VectorWithOffsetType;
				psi_[sectorIndex][excitedIndex]->read(io, label + ttos(excitedIndex));
			}
		}
	}

	const SizeType tvs() const { return targetVectors_.size(); }

	VectorWithOffsetType& targetVectorsNonConst(SizeType i) // <--- FIXME
	{
		assert(i < targetVectors_.size());
		assert(targetVectors_[i]);
		return *targetVectors_[i];
	}

	const VectorWithOffsetType& targetVectors(SizeType i) const
	{
		assert(i < targetVectors_.size());
		assert(targetVectors_[i]);
		return *targetVectors_[i];
	}

	SizeType createPvector(const VectorWithOffsetType& vv)
	{
		const SizeType n = targetVectors_.size();
		VectorWithOffsetType* v = new VectorWithOffsetType();
		targetVectors_.push_back(v);
		*targetVectors_[n] = vv;
		return n;
	}

	void destroyPvector(SizeType ind)
	{
		if (ind >= targetVectors_.size())
			err("AOE: destroyPvector\n");

		targetVectors_[ind]->clear();
	}

	void trimVectors()
	{
		const SizeType n = targetVectors_.size();
		SizeType x = 0;
		for (SizeType ii = 0; ii < n; ++ii) {
			const SizeType i = n - ii - 1;
			if (targetVectors_[i]->size() == 0) {
				++x;
			} else {
				break;
			}
		}

		if (x == 0) return;
		targetVectorsResize(n - x);
	}

	void targetVectorsResize(SizeType x)
	{
		const SizeType n = targetVectors_.size();
		if (n == x) return;

		if (n < x) {
			const SizeType r = x - n;
			for (SizeType i = 0; i < r; ++i) {
				VectorWithOffsetType* v = new VectorWithOffsetType();
				targetVectors_.push_back(v);
			}
		} else {
			assert(n > x);
			const SizeType r = n - x;
			for (SizeType i = 0; i < r; ++i) {
				const SizeType j = n - i - 1;
				delete targetVectors_[j];
				targetVectors_[j] = nullptr;
			}

			targetVectors_.resize(x);
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Number of target vectors is now "<<targetVectors_.size();
		progress_.printline(msgg, std::cout);
	}

	void initTimeVectors(const TargetParamsType& tstStruct,
	                     InputValidatorType& ioIn)
	{
		delete timeVectorsBase_;
		timeVectorsBase_ = nullptr;
		const LeftRightSuperType& lrs = targetHelper_.lrs();
		const ModelType& model = targetHelper_.model();
		const WaveFunctionTransfType& wft = targetHelper_.wft();

		PsimagLite::String s (__FILE__);
		s += " Unknown algorithm\n";

		switch (tstStruct.algorithm()) {
		case TargetParamsType::AlgorithmEnum::KRYLOV:
			timeVectorsBase_ = new TimeVectorsKrylovType(tstStruct,
			                                             targetVectors_,
			                                             model,
			                                             wft,
			                                             lrs,
			                                             ioIn);
			break;
		case TargetParamsType::AlgorithmEnum::CHEBYSHEV:
			timeVectorsBase_ = new TimeVectorsChebyshevType(tstStruct,
			                                                targetVectors_,
			                                                model,
			                                                wft,
			                                                lrs,
			                                                ioIn);
			break;
		case TargetParamsType::AlgorithmEnum::RUNGE_KUTTA:
			timeVectorsBase_ = new TimeVectorsRungeKuttaType(tstStruct,
			                                                 targetVectors_,
			                                                 model,
			                                                 wft,
			                                                 lrs);
			break;
		case TargetParamsType::AlgorithmEnum::SUZUKI_TROTTER:
			timeVectorsBase_ = new TimeVectorsSuzukiTrotterType(tstStruct,
			                                                    targetVectors_,
			                                                    model,
			                                                    wft,
			                                                    lrs);
			break;
		default:
			throw PsimagLite::RuntimeError(s.c_str());
		}
	}

	const TimeVectorsBaseType& timeVectors() const
	{
		assert(timeVectorsBase_);
		return *timeVectorsBase_;
	}

	void loadEnergy(PsimagLite::IoSelector::In& io,
	                PsimagLite::String label)
	{
		SizeType nsectors = 0;
		io.read(nsectors, label + "/Size");
		if (nsectors == 0)
			err("FATAL: " + label + "/Size=0");

		SizeType nexcited = 0;
		try {
			io.read(nexcited, label + "/0/Size");
		} catch (...) {
			loadEnergyLegacy(io, label);
			return;
		}

		if (nexcited == 0)
			err("FATAL: " + label + "/0/Size=0");

		io.read(E0_, label + "/0/0");
	}

	void multiplyTimeVector(SizeType i,RealType factor)
	{
		*targetVectors_[i] = factor*(*targetVectors_[i]);
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
		if (timeVectorsBase_->isBase())
			err("timeVectorsBase_ ptr not setup!?\n");

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

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"I'm applying a local operator now";
		progress_.printline(msgg, std::cout);

		const SizeType splitSize = targetHelper_.model().hilbertSize(site);

		typename PsimagLite::Vector<bool>::Type signs;
		targetHelper_.model().findOddElectronsOfOneSite(signs, site);
		FermionSign fs(targetHelper_.lrs().left(), signs);
		applyOpLocal_(phiNew,
		              phiOld,
		              tstStruct.aOperators()[indexOfOperator],
		              fs,
		              splitSize,
		              systemOrEnviron,
		              corner);

		RealType norma = norm(phiNew);
		if (norma<1e-6) {
			PsimagLite::OstringStream msgg2(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
			msg2<<"Norm of phi is zero\n";
			progress_.printline(msgg2, std::cout);
		}
	}

	void multiSitePush(DmrgSerializerType const* ds) const
	{
		const VectorWithOffsetType& psi00 = ensureOnlyOnePsi("multiSitePush");
		multiSiteExprHelper_.push(ds, psi00);
	}

	void wftSome(SizeType site, SizeType begin, SizeType end)
	{
		wftHelper_.wftSome(targetVectors_, site, begin, end);
	}

	void wftOneVector(VectorWithOffsetType& phiNew,
	                  const VectorWithOffsetType& src,
	                  SizeType site) const
	{
		wftHelper_.wftOneVector(phiNew, src, site);
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

	void timeHasAdvanced()
	{
		timeVectorsBase_->timeHasAdvanced();
	}

	const ModelType& model() const { return targetHelper_.model(); }

private:

	// legacy reading, use only as fallback
	void loadEnergyLegacy(PsimagLite::IoSelector::In& io,
	                      PsimagLite::String label)
	{
		SizeType total = 0;
		io.read(total, label + "/Size");

		for (SizeType i = 0; i < total; ++i)
			io.read(E0_, label + "/" + ttos(i));
	}

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
		const bool advanceOnlyAtBorder = !options.isSet("advanceUnrestricted");

		SizeType sites = targetHelper_.model().superGeometry().numberOfSites();
		bool weAreAtBorder = (site < 1 || site >= sites-1);
		bool dontAdvance = (advanceOnlyAtBorder & !weAreAtBorder);

		if (advanceEach > 0 && timesWithoutAdvancement_ >= advanceEach && !dontAdvance) {
			stage_[i] = StageEnum::WFT_ADVANCE;
			if (i == lastI) {
				timeVectorsBase_->advanceCurrentTimeStep();
				timesWithoutAdvancement_ = 1;
				timeVectorsBase_->timeHasAdvanced();
			}
		} else {
			if (i == lastI &&
			        stage_[i] == StageEnum::WFT_NOADVANCE &&
			        firstSeeLeftCorner_)
				++timesWithoutAdvancement_;
		}

		if (!firstSeeLeftCorner_ &&
		        i==lastI &&
		        stage_[i] == StageEnum::WFT_NOADVANCE &&
		        site==1)
			firstSeeLeftCorner_ = true;

		const RealType time = timeVectorsBase_->time();
		PsimagLite::OstringStream msgg2(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
		msg2<<"Steps without advance: "<<timesWithoutAdvancement_;
		msg2<<" site="<<site<<" currenTime="<<time;
		if (timesWithoutAdvancement_ > 0) progress_.printline(msgg2, std::cout);

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Evolving, stage="<<stageToString(i);
		msg<<" site="<<site<<" loopNumber="<<loopNumber;
		msg<<" Eg="<<Eg;
		progress_.printline(msgg, std::cout);

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

			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"I'm applying a local operator now";
			progress_.printline(msgg, std::cout);
			typename PsimagLite::Vector<bool>::Type signs;
			targetHelper_.model().findOddElectronsOfOneSite(signs, site);

			const SizeType splitSize = targetHelper_.model().hilbertSize(site);
			FermionSign fs(targetHelper_.lrs().left(), signs);
			applyOpLocal_(phiNew,
			              phiOld,
			              tstStruct.aOperators()[i],
			              fs,
			              splitSize,
			              systemOrEnviron,corner);
			RealType norma = norm(phiNew);

			if (norma<1e-6) {
				PsimagLite::OstringStream msgg2(std::cout.precision());
				PsimagLite::OstringStream::OstringStreamType& msg2 = msgg2();
				msg2<<"Norm of phi is zero\n";
				progress_.printline(msgg2, std::cout);
			}
		} else if (stage_[i] == StageEnum::WFT_NOADVANCE ||
		           stage_[i] == StageEnum::WFT_ADVANCE) {
			const SizeType advanceEach = tstStruct.advanceEach();
			SizeType advance = indexNoAdvance_;

			if (advanceEach > 0 &&
			        stage_[i] == StageEnum::WFT_ADVANCE &&
			        !timeVectorsBase_->isBase()) {
				SizeType timeSteps = tstStruct.times().size();
				advance = (timeSteps > 0) ? timeSteps - 1 : 0;
			}

			if (targetVectors_.size() <= advance) {
				PsimagLite::String s(__FILE__);
				s += ": TargetVectors.size()" + ttos(targetVectors_.size());
				s += " but advance=" + ttos(advance) + "\n";
				throw PsimagLite::RuntimeError(s);
			}

			assert(targetVectors_[advance]);
			const VectorWithOffsetType& src = *targetVectors_[advance];

			if (src.size() == 0) {
				PsimagLite::String s(__FILE__);
				s += ": TargetVectors[" + ttos(advance) + "].size()==0\n";
				throw PsimagLite::RuntimeError(s);
			}

			if  (site == 0 || site == numberOfSites -1) {
				// don't wft since we did it before
				assert(advance < targetVectors_.size());
				if (!timeVectorsBase_->isBase() || advanceEach == 0)
					phiNew = src;
				return;
			}

			wftHelper_.wftOneVector(phiNew, src, site);
		} else {
			throw PsimagLite::RuntimeError("computePhi\n");
		}
	}

	void clearPsi()
	{
		const SizeType sectors = psi_.size();
		for (SizeType i = 0; i < sectors; ++i) {
			const SizeType nexcited = psi_[i].size();
			for (SizeType j = 0; j < nexcited; ++j) {
				delete psi_[i][j];
				psi_[i][j] = nullptr;
			}
		}
	}

	ApplyOperatorExpression(const ApplyOperatorExpression&) = delete;

	ApplyOperatorExpression& operator=(const ApplyOperatorExpression&) = delete;

	PsimagLite::ProgressIndicator progress_;
	const TargetHelperType& targetHelper_;
	VectorStageEnumType stage_;
	RealType E0_;
	SizeType indexNoAdvance_;
	ApplyOperatorType applyOpLocal_;
	VectorVectorVectorWithOffsetType psi_;
	VectorVectorWithOffsetType targetVectors_;
	TimeVectorsBaseType* timeVectorsBase_;
	WftHelperType wftHelper_;
	mutable MultiSiteExpressionHelperType multiSiteExprHelper_;
	CorrelationsSkeletonType correlationsSkel_;
	static SizeType timesWithoutAdvancement_;
	static bool firstSeeLeftCorner_;
};

template<typename T1, typename T2, typename T3>
SizeType ApplyOperatorExpression<T1,T2,T3>::timesWithoutAdvancement_ = 0;

template<typename T1, typename T2, typename T3>
bool ApplyOperatorExpression<T1,T2,T3>::firstSeeLeftCorner_ = false;
} // namespace Dmrg

#endif // APPLY_OP_EXPRESSION_H

