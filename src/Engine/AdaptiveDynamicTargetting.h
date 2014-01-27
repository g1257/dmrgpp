/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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

/*! \file AdaptiveDynamicTargetting.h
 *
 * Implements the targetting required by
 * arXiv:1012.5543v1
 *
 */

#ifndef ADAPTIVE_DYN_TARGETTING_H
#define ADAPTIVE_DYN_TARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "DynamicSerializer.h"
#include "AdaptiveDynamicParams.h"
#include "VectorWithOffsets.h"
#include "ParametersForSolver.h"
#include "CommonTargetting.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_,
         typename IoType_>
class AdaptiveDynamicTargetting  {
public:

	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef IoType_ IoType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef AdaptiveDynamicParams<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,MatrixVectorType,VectorType> LanczosSolverType;
	typedef VectorType TargetVectorType;
	typedef typename IoType::In IoInputType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef TargetHelper<ModelType,
	                     TargettingParamsType,
	                     WaveFunctionTransfType> TargetHelperType;
	typedef CommonTargetting<TargetHelperType,
	                         VectorWithOffsetType,
	                         LanczosSolverType> CommonTargettingType;

	enum {DISABLED=CommonTargettingType::DISABLED,
		  OPERATOR=CommonTargettingType::OPERATOR,
		  WFT_NOADVANCE=CommonTargettingType::WFT_NOADVANCE};

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		  EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
		  INFINITE=WaveFunctionTransfType::INFINITE};

	AdaptiveDynamicTargetting(const LeftRightSuperType& lrs,
	                          const ModelType& model,
	                          const TargettingParamsType& tstStruct,
	                          const WaveFunctionTransfType& wft,
	                          const SizeType& quantumSector) // quantumSector ignored here
	    : lrs_(lrs),
	      model_(model),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      lastLanczosVector_(0),
	      dynCounter_(0),
	      progress_("AdaptiveDynamicTargetting"),
	      gsWeight_(1.0),
	      commonTargetting_(lrs,model,tstStruct,wft),
	      done_(false),
	      weightForContinuedFraction_(0)
	{
		if (!wft.isEnabled()) throw PsimagLite::RuntimeError(" DynamicTargetting "
		                                                     "needs an enabled wft\n");
	}

	const ModelType& model() const { return model_; }

	RealType weight(SizeType i) const
	{
		assert(commonTargetting_.allStages(DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		return gsWeight_;
	}

	RealType normSquared(SizeType i) const
	{
		return commonTargetting_.normSquared(i);
	}

	template<typename SomeBasisType>
	void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
	           const SomeBasisType& someBasis)
	{
		commonTargetting_.psi().set(v,someBasis);
	}

	RealType& operator[](SizeType i)
	{
		return commonTargetting_.targetVectors()[i];
	}

	const VectorWithOffsetType& gs() const { return commonTargetting_.psi(); }

	bool includeGroundStage() const
	{
		return (fabs(gsWeight_)>1e-6);
	}

	SizeType size() const
	{
		if (!commonTargetting_.allStages(WFT_NOADVANCE)) return 0;
		return lastLanczosVector_;
	}

	const VectorWithOffsetType& operator()(SizeType i) const
	{
		return commonTargetting_.targetVectors()[i];
	}

	void evolve(RealType Eg,
	            SizeType direction,
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
		SizeType numberOfSites = lrs_.super().block().size();
		if (site>1 && site<numberOfSites-2) return;
		//corner case
		SizeType x = (site==1) ? 0 : numberOfSites-1;
		evolve(Eg,direction,x,loopNumber);
	}

	void initialGuess(VectorWithOffsetType& v,
	                  const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		commonTargetting_.initialGuess(v,block,weight_);
	}

	const LeftRightSuperType& leftRightSuper() const { return lrs_; }

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          IoOutputType& io) const
	{
		assert(block.size()==1);
		SizeType type = tstStruct_.type;
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;

		if (ab_.size()<2) return;
		typename PostProcType::ParametersType params;
		params.Eg = Eg_;
		params.weight = s2*weightForContinuedFraction_;
		params.isign = s;
		PostProcType cf(ab_,reortho_,params);
		PsimagLite::String str = "#TCENTRALSITE=" + ttos(block[0]);
		io.printline(str);
		commonTargetting_.save(block,io,cf,commonTargetting_.targetVectors());

		commonTargetting_.psi().save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		IoInputType io(f);

		try {
			commonTargetting_.setAllStagesTo(WFT_NOADVANCE);

			TimeSerializerType dynS(io,IoInputType::LAST_INSTANCE);
			commonTargetting_.loadTargetVectors(dynS);
			lastLanczosVector_ = commonTargetting_.targetVectors().size()-1;

			//! WARNING: USE OF MAGIC BELOW:
			dynCounter_ = 13; // FIXME: MAYBE SAVE AND LOAD ACTUAL NUMBER HERE
			commonTargetting_.psi().load(io,"PSI");
		} catch (std::exception& e) {
			std::cout<<"WARNING: No special targets found in file "<<f<<"\n";
			commonTargetting_.setAllStagesTo(DISABLED);
			io.rewind();
			int site = 0;
			io.readline(site,"#TCENTRALSITE=",IoType::In::LAST_INSTANCE);
			commonTargetting_.psi().loadOneSector(io,"PSI");
		}
	}

	RealType time() const { return 0; }

	void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
	{
		// nothing to do here
	}

	bool end() const { return false; }

private:

	void evolve(RealType Eg,SizeType direction,SizeType site,SizeType loopNumber)
	{
		Eg_ = Eg;
		VectorWithOffsetType phiNew;
		if (ab_.size()==0)
			commonTargetting_.getPhi(phiNew,Eg,direction,site,loopNumber);

		if (!commonTargetting_.allStages(WFT_NOADVANCE)) {
			commonTargetting_.targetVectors(0) = phiNew;
			return;
		}

		SizeType numberOfSites = lrs_.super().block().size();
		if (site>0 && site<numberOfSites-1) wftAllDynVectors(site);

		if (!done_) calcDynVectors(site,phiNew);
	}

	void wftAllDynVectors(SizeType site)
	{
		for (SizeType i=0;i<=lastLanczosVector_;i++)
			if (i<2) wftOneDynVector(i,site);
	}

	void wftOneDynVector(SizeType i,SizeType site)
	{
		typename PsimagLite::Vector<SizeType>::Type nk(1,model_.hilbertSize(site));

		VectorWithOffsetType result;
		result.populateSectors(lrs_.super());

		// OK, now that we got the partition number right, let's wft:

		// FIXME generalize for su(2)
		wft_.setInitialVector(result,commonTargetting_.targetVectors()[i],lrs_,nk);
		result.collapseSectors();
		commonTargetting_.targetVectors(i) = result;
	}

	void calcDynVectors(SizeType site,const VectorWithOffsetType& phiNew)
	{
		for (SizeType i=0;i<commonTargetting_.targetVectors()[0].sectors();i++) {
			VectorType sv;
			SizeType i0 = commonTargetting_.targetVectors()[0].sector(i);
			commonTargetting_.targetVectors()[0].extract(sv,i0);
			SizeType p = lrs_.super().findPartitionNumber(commonTargetting_.targetVectors()[0].offset(i0));
			if (i==0) {
				if (lastLanczosVector_==0)
					commonTargetting_.targetVectors(1) = commonTargetting_.targetVectors()[0];
			}
			setLanczosVectors(i0,sv,p,site,phiNew);
		}
		setWeights();
		if (lastLanczosVector_==1 && fabs(weightForContinuedFraction_)<1e-6)
			weightForContinuedFraction_ = std::real(commonTargetting_.targetVectors()[0]*commonTargetting_.targetVectors()[0]);
	}

	void setLanczosVectors(
	        SizeType i0,
	        const VectorType& sv,
	        SizeType p,
	        SizeType site,
	        const VectorWithOffsetType& phiNew)
	{
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,threadId);
		typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
		LanczosMatrixType h(&model_,&modelHelper);

		ParametersForSolverType params;
		params.steps = model_.params().lanczosSteps;
		params.tolerance = model_.params().lanczosEps;
		params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

		LanczosSolverType lanczosSolver(h,params);

		RealType a=0,b=0;
		VectorType x(sv.size(),0.0);
		VectorType y = sv;
		if (lastLanczosVector_>0) {
			commonTargetting_.targetVectors()[1].extract(y,i0);
			commonTargetting_.targetVectors()[0].extract(x,i0);
		}

		if (lastLanczosVector_==0) normalize(y);

		lanczosSolver.oneStepDecomposition(x,y,a,b,lastLanczosVector_==0);

		if (lastLanczosVector_<2) lastLanczosVector_++;

		//f0 is wft'd, do nothing
		//f1 is wft'd, do nothing
		RealType norm1 = PsimagLite::norm(x);
		if (norm1<1e-6) {
			lanczosSolver.push(ab_,a,b);
			h.matrixVectorProduct(x,y);
			a = std::real(x*y);
			lanczosSolver.push(ab_,a,b);
			done_=true;
			return;
		}
		dynCounter_++;
		if (lastLanczosVector_>1 && (dynCounter_%tstStruct_.advanceEach) != 0)
			return;
		commonTargetting_.targetVectors(0).setDataInSector(x,i0);
		commonTargetting_.targetVectors(1).setDataInSector(y,i0);
		if ((dynCounter_%tstStruct_.advanceEach) != 0) return;
		if (ab_.size()>0) {
			lanczosSolver.push(ab_,a,b);
			return;
		}
		// first push:
		VectorType xx(sv.size(),0.0);
		VectorType yy;
		phiNew.extract(yy,i0);
		normalize(yy);
		RealType a1=0,b1=0;
		lanczosSolver.oneStepDecomposition(xx,yy,a1,b1,lastLanczosVector_==0);
		lanczosSolver.push(ab_,a1,b1);
		if (tstStruct_.advanceEach<=1) return;
		lanczosSolver.push(ab_,a,b);
	}

	void setWeights()
	{
		RealType sum  = 0;
		weight_.resize(commonTargetting_.targetVectors().size());
		for (SizeType r=0;r<weight_.size();r++) {
			weight_[r] = (r>lastLanczosVector_) ? 0 : 1;
			sum += weight_[r];
		}

		gsWeight_ = commonTargetting_.setGsWeight(0.0);
		for (SizeType r=0;r<weight_.size();r++) weight_[r] *= (1-gsWeight_)/sum;
	}

	RealType dynWeightOf(VectorType& v,const VectorType& w) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<v.size();i++) {
			if (i>=w.size()) continue;
			RealType tmp = std::real(v[i]*w[i]);
			sum += tmp*tmp;
		}
		return sum;
	}

	void normalize(VectorType& v) const
	{
		RealType x = PsimagLite::norm(v);
		if (fabs(x)<1e-6) return;
		v /= x;
	}

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	SizeType lastLanczosVector_;
	SizeType dynCounter_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	CommonTargettingType commonTargetting_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	bool done_;
	RealType Eg_;
	RealType weightForContinuedFraction_;
	TridiagonalMatrixType ab_;
	DenseMatrixRealType reortho_;
}; // class DynamicTargetting

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType,
         typename IoType_>
std::ostream& operator<<(std::ostream& os,
                         const AdaptiveDynamicTargetting<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType,IoType_>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // ADAPTIVE_DYN_TARGETTING_H

