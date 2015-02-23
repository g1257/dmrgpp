/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

#ifndef DMRG_TARGETING_ANCILLA_H
#define DMRG_TARGETING_ANCILLA_H

#include <iostream>
#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TimeSerializer.h"
#include "TargetParamsAncilla.h"
#include "ProgramGlobals.h"
#include "ParametersForSolver.h"
#include "ParallelWft.h"
#include "TimeVectorsKrylov.h"
#include "TimeVectorsRungeKutta.h"
#include "TimeVectorsSuzukiTrotter.h"
#include "TargetingBase.h"
#include "BlockMatrix.h"
#include "MettsStochastics.h"
#include "PackIndices.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_>
class TargetingAncilla : public TargetingBase<LanczosSolverTemplate,
                                               MatrixVectorType_,
                                               WaveFunctionTransfType_> {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	struct AncillaPrev {

		AncillaPrev() : fixed(0),permutationInverse(0)
		{}

		SizeType fixed;
		VectorSizeType permutationInverse;
	};

	enum {BORDER_NEITHER, BORDER_LEFT, BORDER_RIGHT};

	enum {INFINITE=WaveFunctionTransfType_::INFINITE};

	const static SizeType SYSTEM = ProgramGlobals::SYSTEM;

public:

	typedef TargetingBase<LanczosSolverTemplate,
	                      MatrixVectorType_,
	                      WaveFunctionTransfType_> BaseType;
	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType TargetVectorType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              MatrixVectorType,
	                              TargetVectorType> LanczosSolverType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorType;
	typedef PsimagLite::Matrix<ComplexType> ComplexMatrixType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef TargetParamsAncilla<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef BlockMatrix<ComplexMatrixType> ComplexBlockMatrixType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename BasisWithOperatorsType::BasisDataType BasisDataType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::RandomForTests<RealType> RngType;
	typedef MettsStochastics<ModelType,RngType> MettsStochasticsType;
	typedef typename MettsStochasticsType::PairType PairType;

	enum {DISABLED,OPERATOR,WFT_NOADVANCE,WFT_ADVANCE};

	static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
	static SizeType const SUM = TargettingParamsType::SUM;

	TargetingAncilla(const LeftRightSuperType& lrs,
	                  const ModelType& model,
	                  const TargettingParamsType& tstStruct,
	                  const WaveFunctionTransfType& wft,
	                  const SizeType&,
	                  InputValidatorType& ioIn)
	    : BaseType(lrs,model,tstStruct,wft,tstStruct.timeSteps(),0),
	      model_(model),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      progress_("TargetingAncilla"),
	      times_(tstStruct_.timeSteps()),
	      weight_(tstStruct_.timeSteps()),
	      mettsStochastics_(model,tstStruct.rngSeed())
	{
		if (!wft.isEnabled()) {
			PsimagLite::String msg("TargetingAncilla needs an enabled wft\n");
			throw PsimagLite::RuntimeError(msg);
		}

		if (tstStruct_.sites() == 0) {
			PsimagLite::String msg("TargetingAncilla needs at least one TSPSite\n");
			throw PsimagLite::RuntimeError(msg);
		}

		RealType tau =tstStruct_.tau();
		RealType sum = 0;
		SizeType n = times_.size();
		gsWeight_ = (tstStruct_.concatenation() == SUM) ? 0.1 : 0.0;
		gsWeight_ = this->common().setGsWeight(gsWeight_);

		RealType factor = (n+4.0)/(n+2.0);
		factor *= (1.0 - gsWeight_);
		for (SizeType i=0;i<n;i++) {
			times_[i] = i*tau/(n-1);
			weight_[i] = factor/(n+4);
			sum += weight_[i];
		}
		sum -= weight_[0];
		sum -= weight_[n-1];
		weight_[0] = weight_[n-1] = 2*factor/(n+4);
		sum += weight_[n-1];
		sum += weight_[0];

		gsWeight_=1.0-sum;
		sum += gsWeight_;
		assert(fabs(sum-1.0)<1e-5);

		this->common().initTimeVectors(times_,ioIn);
	}

	RealType weight(SizeType i) const
	{
		assert(!this->common().allStages(DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return 0;
	}

	bool includeGroundStage() const
	{
		return false;
	}

	void evolve(RealType Eg,
	            SizeType direction,
	            const BlockType& block1,
	            const BlockType& block2,
	            SizeType loopNumber)
	{
		if (direction==INFINITE) {
			updateStochastics(block1,block2);
			getNewPures(block1,block2);
			return;
		}

		VectorWithOffsetType phiNew;
		this->common().getPhi(phiNew,Eg,direction,block1[0],loopNumber);

		PairType startEnd(0,times_.size());
		bool allOperatorsApplied = (this->common().noStageIs(DISABLED) &&
		                            this->common().noStageIs(OPERATOR));

		this->common().calcTimeVectors(startEnd,
		                                  Eg,
		                                  phiNew,
		                                  direction,
		                                  allOperatorsApplied,
		                                  block1);

		cocoon(direction,block1); // in-situ
		printEnergies(); // in-situ
		printNormsAndWeights();
	}

	void load(const PsimagLite::String& f)
	{
		this->common().template load<TimeSerializerType>(f);
	}

	void print(std::ostream& os) const
	{
		os<<"TSTWeightsTimeVectors=";
		for (SizeType i=0;i<weight_.size();i++)
			os<<weight_[i]<<" ";
		os<<"\n";
		os<<"TSTWeightGroundState="<<gsWeight_<<"\n";
	}

	template<typename IoOutputType>
	void save(const VectorSizeType& block,IoOutputType& io) const
	{
		PsimagLite::OstringStream msg;
		msg<<"Saving state...";
		progress_.printline(msg,std::cout);

		SizeType marker = 0;
		if (this->common().noStageIs(DISABLED)) marker = 1;

		TimeSerializerType ts(this->common().currentTime(),
		                      block[0],
		        this->common().targetVectors(),
		        marker);
		ts.save(io);
		this->common().psi().save(io,"PSI");
	}

	void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
	{

		BlockType X = basisWithOps.block();
		if (X.size()!=1) return;
		assert(X[0]==0 || X[0]==this->leftRightSuper().super().block().size()-1);
		typename PsimagLite::Vector<OperatorType>::Type creationMatrix;
		SparseMatrixType hmatrix;
		BasisDataType q;
		this->model().setNaturalBasis(creationMatrix,
		                              hmatrix,
		                              q,
		                              X,
		                              this->common().currentTime());
		basisWithOps.setVarious(X,hmatrix,q,creationMatrix);
	}

	bool end() const
	{
		return (tstStruct_.maxTime() != 0 &&
		        this->common().currentTime() >= tstStruct_.maxTime());
	}

private:

	void updateStochastics(const VectorSizeType& block1,
	                       const VectorSizeType& block2)
	{
		SizeType linSize = model_.geometry().numberOfSites();
		VectorSizeType tqn(2,0);
		if (model_.params().targetQuantumNumbers.size()>=2) {
			tqn[0] = SizeType(round(model_.params().targetQuantumNumbers[0]*linSize));
			tqn[1] = SizeType(round(model_.params().targetQuantumNumbers[1]*linSize));
		} else {
			tqn[0] = model_.params().electronsUp;
			tqn[1] = model_.params().electronsDown;
		}
		SizeType qn = BasisType::pseudoQuantumNumber(tqn);
		mettsStochastics_.update(qn,block1,block2,tstStruct_.rngSeed());
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

		const SparseMatrixType& transformSystem =  wft_.transform(ProgramGlobals::SYSTEM);
		TargetVectorType newVector1(transformSystem.row(),0);

		VectorSizeType nk1;
		setNk(nk1,block1);
		SizeType alphaFixedVolume = volumeOf(alphaFixed,nk1);

		getNewPure(newVector1,pureVectors_.first,ProgramGlobals::SYSTEM,
		           alphaFixedVolume,this->lrs().left(),transformSystem,block1);
		pureVectors_.first = newVector1;

		const SparseMatrixType& transformEnviron =
		        wft_.transform(ProgramGlobals::ENVIRON);
		TargetVectorType newVector2(transformEnviron.row(),0);

		VectorSizeType nk2;
		setNk(nk2,block2);
		SizeType betaFixedVolume = volumeOf(betaFixed,nk2);
		getNewPure(newVector2,pureVectors_.second,ProgramGlobals::ENVIRON,
		           betaFixedVolume,this->lrs().right(),transformEnviron,block2);
		pureVectors_.second = newVector2;
		setFromInfinite(this->common().targetVectors(0));
		assert(std::norm(targetVectors_[0])>1e-6);

		systemPrev_.fixed = alphaFixedVolume;
		systemPrev_.permutationInverse = this->lrs().left().permutationInverse();
		environPrev_.fixed = betaFixedVolume;
		environPrev_.permutationInverse = this->lrs().right().permutationInverse();
	}

	void getNewPure(TargetVectorType& newVector,
	                TargetVectorType& oldVector,
	                SizeType direction,
	                SizeType alphaFixed,
	                const BasisWithOperatorsType& basis,
	                const SparseMatrixType& transform,
	                const VectorSizeType& block)
	{
		if (oldVector.size()==0)
			setInitialPure(oldVector,block);
		TargetVectorType tmpVector;
		if (transform.row()==0) {
			tmpVector = oldVector;
			assert(PsimagLite::norm(tmpVector)>1e-6);
		} else {
			delayedTransform(tmpVector,oldVector,direction,transform,block);
			assert(PsimagLite::norm(tmpVector)>1e-6);
		}
		SizeType ns = tmpVector.size();
		VectorSizeType nk;
		setNk(nk,block);
		SizeType volumeOfNk = volumeOf(nk);
		SizeType newSize =  (transform.col()==0) ? (ns*ns) :
		                                           transform.col() * volumeOfNk;
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

	void setFromInfinite(VectorWithOffsetType& phi) const
	{
		const LeftRightSuperType& lrs = this->lrs();
		phi.populateSectors(lrs.super());
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i0 = phi.sector(ii);
			TargetVectorType v;
			getFullVector(v,i0,lrs);
			RealType tmpNorm = PsimagLite::norm(v);
			if (fabs(tmpNorm-1.0)<1e-6) {
				SizeType j = lrs.super().qn(lrs.super().partition(i0));
				VectorSizeType qns = BasisType::decodeQuantumNumber(j);
				std::cerr<<"setFromInfinite: qns= ";
				for (SizeType k=0;k<qns.size();k++) std::cerr<<qns[k]<<" ";
				std::cerr<<"\n";
			}
			phi.setDataInSector(v,i0);
		}
		phi.collapseSectors();
		assert(std::norm(phi)>1e-6);
	}

	void setInitialPure(TargetVectorType& oldVector,const VectorSizeType& block)
	{
		int offset = (block[0]==block.size()) ? -block.size() : block.size();
		VectorSizeType blockCorrected = block;
		for (SizeType i=0;i<blockCorrected.size();i++)
			blockCorrected[i] += offset;

		VectorSizeType nk;
		setNk(nk,blockCorrected);
		SizeType volumeOfNk = volumeOf(nk);
		VectorSizeType alphaFixed(nk.size());
		for (SizeType i=0;i<alphaFixed.size();i++)
			alphaFixed[i] = mettsStochastics_.chooseRandomState(blockCorrected[i]);

		PsimagLite::OstringStream msg;
		msg<<"New pures for site ";
		for (SizeType i=0;i<blockCorrected.size();i++)
			msg<<blockCorrected[i]<<" ";
		msg<<" is "<<alphaFixed;
		progress_.printline(msg,std::cerr);

		SizeType volumeOfAlphaFixed = volumeOf(alphaFixed,nk);

		oldVector.resize(volumeOfNk);
		assert(volumeOfAlphaFixed<oldVector.size());
		for (SizeType i=0;i<oldVector.size();i++) {
			oldVector[i] = (i==volumeOfAlphaFixed) ? 1 : 0;
		}
		assert(PsimagLite::norm(oldVector)>1e-6);
	}

	void delayedTransform(TargetVectorType& newVector,
	                      TargetVectorType& oldVector,
	                      SizeType direction,
	                      const SparseMatrixType& transform,
	                      const VectorSizeType& block)
	{
		assert(oldVector.size()==transform.row());

		VectorSizeType nk;
		setNk(nk,block);
		SizeType ne = volumeOf(nk);

		const VectorSizeType& permutationInverse = (direction==SYSTEM)
		        ? systemPrev_.permutationInverse : environPrev_.permutationInverse;
		SizeType nsPrev = permutationInverse.size()/ne;

		newVector.resize(transform.col());
		//newVector = oldVector * transform;
		for (SizeType gamma=0;gamma<newVector.size();gamma++) {
			newVector[gamma] = 0;
			for (SizeType alpha=0;alpha<nsPrev;alpha++) {
				SizeType noPermIndex =  (direction==SYSTEM)
				        ? alpha + systemPrev_.fixed*nsPrev
				        : environPrev_.fixed + alpha*ne;

				SizeType gammaPrime = permutationInverse[noPermIndex];

				assert(gammaPrime<transform.row());
				newVector[gamma] += transform.element(gammaPrime,gamma) *
				        oldVector[gammaPrime];
			}
		}
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

	void setNk(typename PsimagLite::Vector<SizeType>::Type& nk,
	           const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		for (SizeType i=0;i<block.size();i++)
			nk.push_back(mettsStochastics_.model().hilbertSize(block[i]));
	}

	SizeType volumeOf(const typename PsimagLite::Vector<SizeType>::Type& v) const
	{
		assert(v.size()>0);
		SizeType ret = v[0];
		for (SizeType i=1;i<v.size();i++) ret *= v[i];
		return ret;
	}

	SizeType volumeOf(const typename PsimagLite::Vector<SizeType>::Type& alphaFixed,
	                  const typename PsimagLite::Vector<SizeType>::Type& nk) const
	{
		assert(alphaFixed.size()>0);
		assert(alphaFixed.size()==nk.size());
		SizeType sum = alphaFixed[0];
		for (SizeType i=1;i<alphaFixed.size();i++)
			sum += alphaFixed[i]*nk[i-1];
		return sum;
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
		msg2<<"gsNorm="<<std::norm(this->common().psi())<<" norms= ";
		for (SizeType i = 0; i < weight_.size(); i++)
			msg2<<this->common().normSquared(i)<<" ";
		progress_.printline(msg2,std::cout);
	}

	void printEnergies() const
	{
		for (SizeType i=0;i<this->common().targetVectors().size();i++)
			printEnergies(this->common().targetVectors()[i],i);
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
		SizeType p = this->leftRightSuper().super().findPartitionNumber(phi.offset(i0));
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,
		                                                this->leftRightSuper(),
		                                                threadId);
		typename LanczosSolverType::LanczosMatrixType lanczosHelper(&this->model(),
		                                                            &modelHelper);

		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		TargetVectorType x(total);
		lanczosHelper.matrixVectorProduct(x,phi2);
		PsimagLite::OstringStream msg;
		msg<<"Hamiltonian average at time="<<this->common().currentTime();
		msg<<" for target="<<whatTarget;
		msg<<" sector="<<i0<<" <phi(t)|H|phi(t)>="<<(phi2*x);
		msg<<" <phi(t)|phi(t)>="<<(phi2*phi2);
		progress_.printline(msg,std::cout);
	}

	// in situ computation:
	void cocoon(SizeType direction,const BlockType& block) const
	{
		std::cout<<"-------------&*&*&* In-situ measurements start\n";

		if (this->common().noStageIs(DISABLED))
			std::cout<<"ALL OPERATORS HAVE BEEN APPLIED\n";
		else
			std::cout<<"NOT ALL OPERATORS APPLIED YET\n";

		PsimagLite::String modelName = this->model().params().model;

		if (modelName == "HubbardOneBand" ||
		    modelName == "HubbardOneBandExtended" ||
		    modelName == "Immm") {
			this->common().cocoonLegacy(direction,block);
		}

		this->common().cocoon(block,direction);

		std::cout<<"-------------&*&*&* In-situ measurements end\n";
	}

	const ModelType& model_;
	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	typename PsimagLite::Vector<RealType>::Type times_,weight_;
	RealType gsWeight_;
	MettsStochasticsType mettsStochastics_;
	AncillaPrev systemPrev_;
	AncillaPrev environPrev_;
	std::pair<TargetVectorType,TargetVectorType> pureVectors_;
};     //class TargetingAncilla

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingAncilla<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType>& tst)
{
	tst.print(os);
	return os;
}
} // namespace Dmrg

#endif

