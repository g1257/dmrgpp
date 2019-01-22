/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

#ifndef DMRG_TARGETING_METTS_H
#define DMRG_TARGETING_METTS_H
#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "MettsSerializer.h"
#include "MettsParams.h"
#include "MettsStochastics.h"
#include <cassert>
#include "MettsCollapse.h"
#include "VectorWithOffset.h"
#include "ParametersForSolver.h"
#include "RandomForTests.h"
#include "TimeSerializer.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include "CrsMatrix.h"
#include "TargetingBase.h"
#include "Io/IoSelector.h"

namespace Dmrg {

template<typename LanczosSolverType_, typename VectorWithOffsetType_>
class TargetingMetts  : public TargetingBase<LanczosSolverType_,VectorWithOffsetType_> {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

	struct MettsPrev {

		MettsPrev() : fixed(0),permutationInverse(0)
		{}

		SizeType fixed;
		VectorSizeType permutationInverse;
	};

public:

	typedef LanczosSolverType_ LanczosSolverType;
	typedef TargetingBase<LanczosSolverType,VectorWithOffsetType_> BaseType;
	typedef typename BaseType::TargetingCommonType TargetingCommonType;
	typedef typename BaseType::MatrixVectorType MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename BaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType TargetVectorType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef MettsParams<ModelType> TargetParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef ApplyOperatorLocal<LeftRightSuperType,VectorWithOffsetType> ApplyOperatorType;
	typedef typename ApplyOperatorType::BorderEnum BorderEnumType;
	typedef MettsSerializer<VectorWithOffsetType> MettsSerializerType;
	typedef typename PsimagLite::RandomForTests<RealType> RngType;
	typedef MettsStochastics<ModelType,RngType> MettsStochasticsType;
	typedef typename MettsStochasticsType::PairType PairType;
	typedef MettsCollapse<VectorWithOffsetType,
	MettsStochasticsType,
	TargetParamsType> MettsCollapseType;
	typedef typename MettsCollapseType::PackIndicesType PackIndicesType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef TimeVectorsBase<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsBaseType;
	typedef TimeVectorsKrylov<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsKrylovType;
	typedef TimeVectorsRungeKutta<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsRungeKuttaType;
	typedef TimeVectorsSuzukiTrotter<TargetParamsType,ModelType,WaveFunctionTransfType,
	LanczosSolverType,VectorWithOffsetType> TimeVectorsSuzukiTrotterType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename BasisType::QnType QnType;
	typedef typename PsimagLite::Vector<BlockDiagonalMatrixType*>::Type
	VectorBlockDiagonalMatrixType;
	typedef typename TargetingCommonType::ApplyOperatorExpressionType::StageEnum StageEnumType;

	const static SizeType SYSTEM = ProgramGlobals::SYSTEM;

	TargetingMetts(const LeftRightSuperType& lrs,
	               const ModelType& model,
	               const WaveFunctionTransfType& wft,
	               const QnType& quantumSector,
	               InputValidatorType& ioIn)
	    : BaseType(lrs,model,wft,0),
	      model_(model),
	      lrs_(lrs),
	      mettsStruct_(ioIn,model),
	      wft_(wft),
	      quantumSector_(quantumSector),
	      progress_("TargetingMetts"),
	      mettsStochastics_(model,mettsStruct_.rngSeed,mettsStruct_.pure),
	      mettsCollapse_(mettsStochastics_,lrs,mettsStruct_),
	      prevDirection_(ProgramGlobals::INFINITE),
	      systemPrev_(),
	      environPrev_()
	{
		if (!wft.isEnabled()) err(" TargetingMetts needs an enabled wft\n");

		RealType tau = mettsStruct_.tau()/(mettsStruct_.timeSteps()-1);
		SizeType n1 = mettsStruct_.timeSteps();
		SizeType n = mettsStruct_.timeSteps() + 1;

		betas_.resize(n1);
		weight_.resize(n);

		gsWeight_= 0.0;
		RealType factor = (1.0 - gsWeight_)/(n1+4);
		RealType sum = setOneInterval(factor,PairType(0,n1),tau);
		weight_[n-1] = 2*factor;
		sum += weight_[n-1];

		sum += gsWeight_;
		assert(fabs(sum-1.0)<1e-5);

		this->common().aoe().initTimeVectors(mettsStruct_, betas_, ioIn);
	}

	~TargetingMetts()
	{
		SizeType n = garbage_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	SizeType sites() const { return mettsStruct_.sites(); }

	SizeType targets() const { return mettsStruct_.timeSteps() + 1; }

	RealType weight(SizeType i) const
	{
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeight_;
	}

	bool includeGroundStage() const {return (fabs(gsWeight_)>1e-6); }

	SizeType size() const
	{
		if (this->common().aoe().allStages(StageEnumType::DISABLED))
			return 1;
		SizeType n = this->common().aoe().targetVectors().size();
		if (this->common().aoe().targetVectors()[n-1].size()==0) n--;
		return n;
	}

	void evolve(RealType Eg,
	            ProgramGlobals::DirectionEnum direction,
	            const BlockType& block1,
	            const BlockType& block2,
	            SizeType loopNumber)
	{
		VectorSizeType sites;
		if (direction == ProgramGlobals::INFINITE)
			utils::blockUnion(sites,block1,block2);
		else sites = block1;

		SizeType n1 = mettsStruct_.timeSteps();

		if (direction == ProgramGlobals::INFINITE) {
			updateStochastics(block1,block2);
			getNewPures(block1,block2);
			return;
		}

		SizeType max = n1;

		if (this->common().aoe().noStageIs(StageEnumType::DISABLED)) {
			max = 1;
			if (this->common().aoe().allStages(StageEnumType::WFT_ADVANCE))
				this->common().aoe().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
		}

		// Advance or wft each target vector for beta/2
		for (SizeType i=0;i<max;i++) {
			evolve(i,0,n1-1,Eg,direction,sites,loopNumber);
		}

		// compute imag. time evolution:
		calcTimeVectors(PairType(0,n1),Eg,direction,block1);

		// Advance or wft  collapsed vector
		if (this->common().aoe().targetVectors()[n1].size()>0)
			evolve(n1,n1,n1-1,Eg,direction,sites,loopNumber);

		for (SizeType i=0;i<this->common().aoe().targetVectors().size();i++)
			assert(this->common().aoe().targetVectors()[i].size()==0 ||
			       this->common().aoe().targetVectors()[i].size()==
			       lrs_.super().permutationVector().size());

		this->common().cocoon(block1, direction);

		printEnergies(); // in-situ

		PsimagLite::String options = this->model().params().options;
		bool normalizeTimeVectors = true;
		if (options.find("neverNormalizeVectors") != std::string::npos)
			normalizeTimeVectors = false;

		if (normalizeTimeVectors)
			this->common().normalizeTimeVectors();

		printNormsAndWeights();

		if (this->common().aoe().noStageIs(StageEnumType::COLLAPSE)) return;

		// collapse
		bool hasCollapsed = mettsCollapse_(this->common().aoe().targetVectors(n1),
		                                   this->common().aoe().targetVectors()[n1-1],
		        sites,
		        direction);
		if (hasCollapsed) {
			PsimagLite::OstringStream msg;
			msg<<"Has Collapsed";
			progress_.printline(msg,std::cout);
		}
	}

	void read(typename TargetingCommonType::IoInputType& io, PsimagLite::String prefix)
	{
		this->common().template readGSandNGSTs<TimeSerializerType>(io, prefix);
	}

	void write(const VectorSizeType& block,
	           PsimagLite::IoSelector::Out& io,
	           PsimagLite::String prefix) const
	{
		this->common().write(io, block, prefix);

		VectorVectorWithOffsetType& tv = const_cast<VectorVectorWithOffsetType&>
		        (this->common().aoe().targetVectors());
		if (mettsStruct_.beta > this->common().aoe().currentTime()) {
			for (SizeType i = 0; i < this->common().aoe().targetVectors().size(); ++i)
				tv[i].clear();
		}

		this->common().writeNGSTs(io, block, prefix);
	}

	void updateOnSiteForCorners(BasisWithOperatorsType&) const
	{
		// nothing to do here
	}

	bool end() const { return false; }

private:

	void evolve(SizeType index,
	            SizeType start,
	            SizeType indexAdvance,
	            RealType Eg,
	            SizeType direction,
	            const VectorSizeType& block,
	            SizeType loopNumber)
	{
		if (index==0 && start==0)
			advanceCounterAndComputeStage(block);

		PsimagLite::OstringStream msg;
		msg<<"Evolving, stage="<<getStage()<<" loopNumber="<<loopNumber;
		msg<<" Eg="<<Eg;
		progress_.printline(msg,std::cout);
		advanceOrWft(index,indexAdvance,direction,block);
	}

	void calcTimeVectors(const PairType& startEnd,
	                     RealType Eg,
	                     ProgramGlobals::DirectionEnum systemOrEnviron,
	                     const VectorSizeType& block)
	{
		const VectorWithOffsetType& phi = this->common().aoe().targetVectors()[startEnd.first];
		PsimagLite::OstringStream msg;
		msg<<" vector number "<<startEnd.first<<" has norm ";
		msg<<norm(phi);
		progress_.printline(msg,std::cout);
		if (norm(phi)<1e-6)
			setFromInfinite(this->common().aoe().targetVectors(startEnd.first),lrs_);
		bool allOperatorsApplied = (this->common().aoe().noStageIs(StageEnumType::DISABLED));
		this->common().aoe().calcTimeVectors(startEnd,
		                               Eg,
		                               phi,
		                               systemOrEnviron,
		                               allOperatorsApplied,
		                               block,
		                               mettsStruct_);
		this->common().normalizeTimeVectors(startEnd.first+1,startEnd.second);
	}

	void advanceCounterAndComputeStage(const VectorSizeType& block)
	{
		static SizeType timesWithoutAdvancement = 0;

		if (this->common().aoe().noStageIs(StageEnumType::COLLAPSE))
			this->common().aoe().setAllStagesTo(StageEnumType::WFT_NOADVANCE);

		if (this->common().aoe().allStages(StageEnumType::COLLAPSE)) {
			if (!allSitesCollapsed()) {
				if (sitesCollapsed_.size()>2*model_.geometry().numberOfSites())
					throw PsimagLite::RuntimeError("advanceCounterAndComputeStage\n");
				printAdvancement(timesWithoutAdvancement);
				return;
			}

			sitesCollapsed_.clear();
			this->common().aoe().setAllStagesTo(StageEnumType::WFT_NOADVANCE);
			timesWithoutAdvancement = 0;
			this->common().aoe().setTime(0);
			PsimagLite::OstringStream msg;
			SizeType n1 = mettsStruct_.timeSteps();
			RealType x = norm(this->common().aoe().targetVectors()[n1]);
			msg<<"Changing direction, setting collapsed with norm="<<x;
			progress_.printline(msg,std::cout);
			for (SizeType i=0;i<n1;i++)
				this->common().aoe().targetVectors(i) = this->common().aoe().targetVectors()[n1];
			this->common().aoe().timeHasAdvanced();
			printAdvancement(timesWithoutAdvancement);
			return;
		}

		if (timesWithoutAdvancement < mettsStruct_.advanceEach()) {
			timesWithoutAdvancement++;
			printAdvancement(timesWithoutAdvancement);
			return;
		}

		if (this->common().aoe().noStageIs(StageEnumType::COLLAPSE) &&
		        this->common().aoe().currentTime() < mettsStruct_.beta) {
			this->common().aoe().setAllStagesTo(StageEnumType::WFT_ADVANCE);
			RealType tmp = this->common().aoe().currentTime() + mettsStruct_.tau();
			this->common().aoe().setTime(tmp);
			timesWithoutAdvancement = 0;
			printAdvancement(timesWithoutAdvancement);
			return;
		}

		if (this->common().aoe().noStageIs(StageEnumType::COLLAPSE) &&
		        this->common().aoe().currentTime() >= mettsStruct_.beta &&
		        block[0]!=block.size()) {
			printAdvancement(timesWithoutAdvancement);
			return;
		}

		if (this->common().aoe().noStageIs(StageEnumType::COLLAPSE) &&
		        this->common().aoe().currentTime() >= mettsStruct_.beta) {
			this->common().aoe().setAllStagesTo(StageEnumType::COLLAPSE);
			sitesCollapsed_.clear();
			SizeType n1 = mettsStruct_.timeSteps();
			this->common().aoe().targetVectors(n1).clear();
			timesWithoutAdvancement = 0;
			printAdvancement(timesWithoutAdvancement);
			return;
		}
	}

	void printAdvancement(SizeType timesWithoutAdvancement) const
	{
		PsimagLite::OstringStream msg2;
		msg2<<"Steps without advance: "<<timesWithoutAdvancement;
		if (timesWithoutAdvancement > 0)
			progress_.printline(msg2,std::cout);
	}

	void advanceOrWft(SizeType index,
	                  SizeType indexAdvance,
	                  SizeType,
	                  const VectorSizeType& block)
	{
		if (this->common().aoe().targetVectors()[index].size()==0) return;
		assert(norm(this->common().aoe().targetVectors()[index])>1e-6);
		VectorSizeType nk;
		mettsCollapse_.setNk(nk,block);

		if (this->common().aoe().allStages(StageEnumType::WFT_NOADVANCE) ||
		        this->common().aoe().allStages(StageEnumType::WFT_ADVANCE) ||
		        this->common().aoe().allStages(StageEnumType::COLLAPSE)) {
			SizeType advance = index;
			if (this->common().aoe().allStages(StageEnumType::WFT_ADVANCE)) {
				advance = indexAdvance;
				this->common().aoe().timeHasAdvanced();
			}
			// don't advance the collapsed vector because we'll recompute
			if (index==weight_.size()-1) advance=index;
			PsimagLite::OstringStream msg;
			msg<<"I'm calling the WFT now";
			progress_.printline(msg,std::cout);

			VectorWithOffsetType phiNew; // same sectors as g.s.
			//phiNew.populateSectors(lrs_.super());
			assert(norm(this->common().aoe().targetVectors()[advance])>1e-6);

			phiNew.populateSectors(lrs_.super());
			// OK, now that we got the partition number right, let's wft:
			wft_.setInitialVector(phiNew,this->common().aoe().targetVectors()[advance],lrs_,nk);
			phiNew.collapseSectors();
			assert(norm(phiNew)>1e-6);
			this->common().aoe().targetVectors(index) = phiNew;
		} else {
			assert(false);
		}
	}

	void updateStochastics(const VectorSizeType& block1,
	                       const VectorSizeType& block2)
	{
		const QnType& qn = model_.targetQuantum().qn;
		mettsStochastics_.update(qn,block1,block2,mettsStruct_.rngSeed);
	}

	SizeType getPartition() const
	{
		SizeType total = lrs_.super().partition()-1;
		for (SizeType i=0;i<total;i++) {
			// Do only one sector unless doing su(2) with j>0, then do all m's
			if (lrs_.super().pseudoQnEqual(i, quantumSector_)) return i;
		}
		throw PsimagLite::RuntimeError("TargetingMetts: getPartition()\n");
	}

	// direction here is INFINITE
	void getNewPures(const VectorSizeType& block1,
	                 const VectorSizeType& block2)
	{
		VectorSizeType alphaFixed(block1.size());
		for (SizeType i=0;i<alphaFixed.size();i++)
			alphaFixed[i] = mettsStochastics_.chooseRandomState(block1[i]);

		VectorSizeType betaFixed(block2.size());
		for (SizeType i=0;i<betaFixed.size();i++)
			betaFixed[i] = mettsStochastics_.chooseRandomState(block2[i]);

		PsimagLite::OstringStream msg;
		msg<<"New pures for ";
		for (SizeType i=0;i<alphaFixed.size();i++)
			msg<<" site="<<block1[i]<<" is "<<alphaFixed[i];
		msg<<" and for ";
		for (SizeType i=0;i<betaFixed.size();i++)
			msg<<" site="<<block2[i]<<" is "<<betaFixed[i];
		progress_.printline(msg,std::cerr);

		const BlockDiagonalMatrixType& transformSystem = getTransform(ProgramGlobals::SYSTEM);

		TargetVectorType newVector1(transformSystem.rows(),0);

		VectorSizeType nk1;
		mettsCollapse_.setNk(nk1,block1);
		SizeType alphaFixedVolume = mettsCollapse_.volumeOf(alphaFixed,nk1);

		getNewPure(newVector1,
		           pureVectors_.first,
		           ProgramGlobals::SYSTEM,
		           alphaFixedVolume,
		           lrs_.left(),
		           transformSystem,
		           block1);
		pureVectors_.first = newVector1;

		const BlockDiagonalMatrixType& transformEnviron = getTransform(ProgramGlobals::ENVIRON);
		TargetVectorType newVector2(transformEnviron.rows(),0);

		VectorSizeType nk2;
		mettsCollapse_.setNk(nk2,block2);
		SizeType betaFixedVolume = mettsCollapse_.volumeOf(betaFixed,nk2);
		getNewPure(newVector2,
		           pureVectors_.second,
		           ProgramGlobals::ENVIRON,
		           betaFixedVolume,
		           lrs_.right(),
		           transformEnviron,
		           block2);
		pureVectors_.second = newVector2;
		setFromInfinite(this->common().aoe().targetVectors(0),lrs_);
		assert(norm(this->common().aoe().targetVectors()[0])>1e-6);

		systemPrev_.fixed = alphaFixedVolume;
		systemPrev_.permutationInverse = lrs_.left().permutationInverse();
		environPrev_.fixed = betaFixedVolume;
		environPrev_.permutationInverse = lrs_.right().permutationInverse();
	}

	void getFullVector(TargetVectorType& v,
	                   SizeType m,
	                   const LeftRightSuperType& lrs) const
	{
		int offset = lrs.super().partition(m);
		int total = lrs.super().partition(m+1) - offset;

		PackIndicesType pack(lrs.left().size());
		v.resize(total);
		assert(PsimagLite::norm(pureVectors_.first)>1e-6);
		assert(PsimagLite::norm(pureVectors_.second)>1e-6);
		for (int i=0;i<total;i++) {
			SizeType alpha,beta;
			pack.unpack(alpha,beta,lrs.super().permutation(i+offset));
			v[i] = pureVectors_.first[alpha] * pureVectors_.second[beta];
		}
	}

	void getNewPure(TargetVectorType& newVector,
	                TargetVectorType& oldVector,
	                SizeType direction,
	                SizeType alphaFixed,
	                const BasisWithOperatorsType& basis,
	                const BlockDiagonalMatrixType& transform,
	                const VectorSizeType& block)
	{
		if (oldVector.size()==0)
			setInitialPure(oldVector,block);
		TargetVectorType tmpVector;
		if (transform.rows()==0) {
			tmpVector = oldVector;
			assert(PsimagLite::norm(tmpVector)>1e-6);
		} else {
			MatrixType transform1;
			transform.toDense(transform1);
			delayedTransform(tmpVector,oldVector,direction,transform1,block);
			assert(PsimagLite::norm(tmpVector)>1e-6);
		}
		SizeType ns = tmpVector.size();
		VectorSizeType nk;
		mettsCollapse_.setNk(nk,block);
		SizeType volumeOfNk = mettsCollapse_.volumeOf(nk);
		SizeType newSize =  (transform.cols()==0) ? (ns*ns) :
		                                            transform.cols() * volumeOfNk;
		newVector.resize(newSize);
		for (SizeType alpha=0;alpha<newVector.size();alpha++)
			newVector[alpha] = 0;

		for (SizeType alpha=0;alpha<ns;alpha++) {
			SizeType gamma = (direction==ProgramGlobals::SYSTEM) ?
			            basis.permutationInverse(alpha + alphaFixed*ns) :
			            basis.permutationInverse(alphaFixed + alpha*volumeOfNk);
			newVector[gamma] = tmpVector[alpha];
		}

		PsimagLite::OstringStream msg2;
		msg2<<"Old size of pure is "<<ns<<" norm="<<PsimagLite::norm(tmpVector);
		progress_.printline(msg2,std::cerr);
		PsimagLite::OstringStream msg;
		msg<<"New size of pure is "<<newSize<<" norm="<<PsimagLite::norm(newVector);
		progress_.printline(msg,std::cerr);
		assert(PsimagLite::norm(newVector)>1e-6);
	}

	void delayedTransform(TargetVectorType& newVector,
	                      TargetVectorType& oldVector,
	                      SizeType direction,
	                      const MatrixType& transform,
	                      const VectorSizeType& block)
	{
		assert(oldVector.size()==transform.rows());

		VectorSizeType nk;
		mettsCollapse_.setNk(nk,block);
		SizeType ne = mettsCollapse_.volumeOf(nk);

		const VectorSizeType& permutationInverse = (direction==SYSTEM)
		        ? systemPrev_.permutationInverse : environPrev_.permutationInverse;
		SizeType nsPrev = permutationInverse.size()/ne;

		newVector.resize(transform.cols());
		//newVector = oldVector * transform;
		for (SizeType gamma=0;gamma<newVector.size();gamma++) {
			newVector[gamma] = 0;
			for (SizeType alpha=0;alpha<nsPrev;alpha++) {
				SizeType noPermIndex =  (direction==SYSTEM)
				        ? alpha + systemPrev_.fixed*nsPrev
				        : environPrev_.fixed + alpha*ne;

				SizeType gammaPrime = permutationInverse[noPermIndex];

				assert(gammaPrime<transform.rows());
				newVector[gamma] += transform(gammaPrime,gamma) *
				        oldVector[gammaPrime];
			}
		}
	}

	void setInitialPure(TargetVectorType& oldVector,const VectorSizeType& block)
	{
		int offset = (block[0]==block.size()) ? -block.size() : block.size();
		VectorSizeType blockCorrected = block;
		for (SizeType i=0;i<blockCorrected.size();i++)
			blockCorrected[i] += offset;

		VectorSizeType nk;
		mettsCollapse_.setNk(nk,blockCorrected);
		SizeType volumeOfNk = mettsCollapse_.volumeOf(nk);
		VectorSizeType alphaFixed(nk.size());
		for (SizeType i=0;i<alphaFixed.size();i++)
			alphaFixed[i] = mettsStochastics_.chooseRandomState(blockCorrected[i]);

		PsimagLite::OstringStream msg;
		msg<<"New pures for site ";
		for (SizeType i=0;i<blockCorrected.size();i++)
			msg<<blockCorrected[i]<<" ";
		msg<<" is "<<alphaFixed;
		progress_.printline(msg,std::cerr);

		SizeType volumeOfAlphaFixed = mettsCollapse_.volumeOf(alphaFixed,nk);

		oldVector.resize(volumeOfNk);
		assert(volumeOfAlphaFixed<oldVector.size());
		for (SizeType i=0;i<oldVector.size();i++) {
			oldVector[i] = (i==volumeOfAlphaFixed) ? 1 : 0;
		}

		assert(PsimagLite::norm(oldVector)>1e-6);
	}

	void setFromInfinite(VectorWithOffsetType& phi,
	                     const LeftRightSuperType& lrs) const
	{
		phi.populateSectors(lrs.super());
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i0 = phi.sector(ii);
			TargetVectorType v;
			getFullVector(v,i0,lrs);
			RealType tmpNorm = PsimagLite::norm(v);
			if (fabs(tmpNorm-1.0)<1e-6) {
				const QnType& j = lrs.super().qnEx(i0);
				std::cerr<<"setFromInfinite: qns= ";
				std::cerr<<j<<"\n";
			}

			phi.setDataInSector(v,i0);
		}
		phi.collapseSectors();
		assert(norm(phi)>1e-6);
	}

	PsimagLite::String getStage() const
	{
		if (this->common().aoe().allStages(StageEnumType::DISABLED))
			return "Disabled";
		if (this->common().aoe().allStages(StageEnumType::COLLAPSE))
			return "Collapsing";
		if (this->common().aoe().allStages(StageEnumType::WFT_ADVANCE))
			return "WFT with time stepping";
		if (this->common().aoe().allStages(StageEnumType::WFT_NOADVANCE))
			return "WFT without time change";

		return "undefined";
	}

	RealType setOneInterval(const RealType& factor,
	                        const PairType& startEnd,
	                        const RealType& tau)
	{
		RealType sum = 0;
		for (SizeType i=startEnd.first;i<startEnd.second;i++) {
			betas_[i] = (i-startEnd.first)*tau;
			weight_[i] = factor;
			sum += weight_[i];
		}

		sum -= weight_[startEnd.first];
		sum -= weight_[startEnd.second-1];
		weight_[startEnd.first] = weight_[startEnd.second-1] = 2*factor;
		sum += weight_[startEnd.second-1];
		sum += weight_[startEnd.first];
		return sum;
	}

	void findElectronsOfOneSite(VectorSizeType& electrons,SizeType site) const
	{
		VectorSizeType block(1,site);
		typename ModelType::HilbertBasisType basis;
		VectorSizeType quantumNumbs;
		model_.setNaturalBasis(basis,quantumNumbs,block);
		model_.findElectrons(electrons,basis,site);
	}

	bool allSitesCollapsed() const
	{
		SizeType n = model_.geometry().numberOfSites();
		for (SizeType i=0;i<n;i++) {
			bool seen = (std::find(sitesCollapsed_.begin(),
			                       sitesCollapsed_.end(),
			                       i) != sitesCollapsed_.end());
			if (!seen) return false;
		}
		return true;
	}

	void printNormsAndWeights() const
	{
		if (this->common().aoe().allStages(StageEnumType::DISABLED))
			return;

		PsimagLite::OstringStream msg;
		msg<<"gsWeight="<<gsWeight_<<" weights= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg<<weight_[i]<<" ";
		progress_.printline(msg,std::cout);

		PsimagLite::OstringStream msg2;
		msg2<<"gsNorm="<<norm(this->common().aoe().psi())<<" norms= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg2<<this->common().normSquared(i)<<" ";
		progress_.printline(msg2,std::cout);
	}

	void printEnergies() const
	{
		for (SizeType i=0;i<this->common().aoe().targetVectors().size();i++)
			printEnergies(this->common().aoe().targetVectors()[i],i);
	}

	void printEnergies(const VectorWithOffsetType& phi,SizeType whatTarget) const
	{
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i = phi.sector(ii);
			printEnergies(phi,whatTarget,i);
		}
	}

	void printEnergies(const VectorWithOffsetType& phi,
	                   SizeType whatTarget,
	                   SizeType i0) const
	{
		SizeType p = this->lrs().super().findPartitionNumber(phi.offset(i0));
		typename ModelType::HamiltonianConnectionType hc(p,
		                                                 BaseType::lrs(),
		                                                 BaseType::model().geometry(),
		                                                 BaseType::ModelType::modelLinks(),
		                                                 this->common().aoe().currentTime(),
		                                                 0);
		typename LanczosSolverType::MatrixType lanczosHelper(BaseType::model(),
		                                                     hc);

		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x,phi2);
		PsimagLite::OstringStream msg;
		msg<<"Hamiltonian average at time="<<this->common().aoe().currentTime();
		msg<<" for target="<<whatTarget;
		ComplexOrRealType numerator = phi2*x;
		ComplexOrRealType den = phi2*phi2;
		ComplexOrRealType division = (PsimagLite::norm(den)<1e-10) ? 0 : numerator/den;
		msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<numerator;
		msg<<" <phi(t)|phi(t)>="<<den<<" "<<division;
		progress_.printline(msg,std::cout);
	}

	const BlockDiagonalMatrixType& getTransform(ProgramGlobals::SysOrEnvEnum sysOrEnv)
	{
		const SizeType stackSize = wft_.size(sysOrEnv);

		if (stackSize > 0) return wft_.getTransform(sysOrEnv);

		BlockDiagonalMatrixType* m = new BlockDiagonalMatrixType;
		garbage_.push_back(m);

		return *m;
	}

	const ModelType& model_;
	const LeftRightSuperType& lrs_;
	TargetParamsType mettsStruct_;
	const WaveFunctionTransfType& wft_;
	const QnType& quantumSector_;
	PsimagLite::ProgressIndicator progress_;
	VectorRealType betas_;
	VectorRealType weight_;
	RealType gsWeight_;
	MettsStochasticsType mettsStochastics_;
	MettsCollapseType mettsCollapse_;
	SizeType prevDirection_;
	MettsPrev systemPrev_;
	MettsPrev environPrev_;
	std::pair<TargetVectorType,TargetVectorType> pureVectors_;
	VectorSizeType sitesCollapsed_;
	VectorBlockDiagonalMatrixType garbage_;
};     //class TargetingMetts
} // namespace Dmrg

#endif //DMRG_TARGETING_METTS_H

