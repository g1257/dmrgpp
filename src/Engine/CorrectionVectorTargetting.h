/*
Copyright (c) 2009-2011, UT-Battelle, LLC
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

/*! \file CorrectionVectorTargetting.h
 *
 * Implements the targetting required by
 * the correction targetting method
 *
 */

#ifndef CORRECTION_VECTOR_TARG_H
#define CORRECTION_VECTOR_TARG_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ApplyOperatorLocal.h"
#include "CorrectionVectorParams.h"
#include "VectorWithOffsets.h"
#include "CorrectionVectorFunction.h"
#include "CommonTargetting.h"
#include "ParametersForSolver.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         template<typename,typename> class InternalProductTemplate,
         template<typename,typename> class WaveFunctionTransfTemplate,
         typename ModelType_,
         typename IoType_,
         template<typename> class VectorWithOffsetTemplate>
class CorrectionVectorTargetting  {

public:

	typedef ModelType_ ModelType;
	typedef IoType_ IoType;
	typedef typename ModelType::RealType RealType;
	typedef std::complex<RealType> ComplexType;
	typedef InternalProductTemplate<RealType,ModelType> InternalProductType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType
	LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef CorrectionVectorParams<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef VectorWithOffsetTemplate<RealType> VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              InternalProductType,
	                              VectorType> LanczosSolverType;
	typedef VectorType TargetVectorType;
	typedef ApplyOperatorLocal<LeftRightSuperType,
	                           VectorWithOffsetType> ApplyOperatorType;
	typedef typename ApplyOperatorType::BorderEnum BorderEnumType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef WaveFunctionTransfTemplate<LeftRightSuperType,
	                                   VectorWithOffsetType> WaveFunctionTransfType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef DynamicSerializer<VectorWithOffsetType,PostProcType> DynamicSerializerType;
	typedef typename LanczosSolverType::LanczosMatrixType
	LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,
	TargettingParamsType>
	CorrectionVectorFunctionType;
	typedef CommonTargetting<ModelType,
	                         TargettingParamsType,
	                         WaveFunctionTransfType,
	                         VectorWithOffsetType,
	                         LanczosSolverType> CommonTargettingType;

	enum {DISABLED,OPERATOR,CONVERGING};

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
	      EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
	      INFINITE=WaveFunctionTransfType::INFINITE};

	static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
	static SizeType const SUM = TargettingParamsType::SUM;

	CorrectionVectorTargetting(const LeftRightSuperType& lrs,
	                           const ModelType& model,
	                           const TargettingParamsType& tstStruct,
	                           const WaveFunctionTransfType& wft,
	                           const SizeType& quantumSector) // quantumSector ignored here
	    : stage_(tstStruct.sites.size(),DISABLED),
	      lrs_(lrs),
	      model_(model),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      progress_("CorrectionVectorTargetting"),
	      applyOpLocal_(lrs),
	      gsWeight_(1.0),
	      targetVectors_(4),
	      commonTargetting_(lrs,model,tstStruct),
	      correctionEnabled_(false)
	{
		if (!wft.isEnabled())
			throw PsimagLite::RuntimeError("CorrectionVectorTargetting needs wft\n");
	}

	const ModelType& model() const { return model_; }

	RealType weight(SizeType i) const
	{
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (!correctionEnabled_) return 1.0;
		return gsWeight_;
	}

	RealType normSquared(SizeType i) const
	{
		return commonTargetting_.normSquared(targetVectors_[i]);
	}

	template<typename SomeBasisType>
	void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
	           const SomeBasisType& someBasis)
	{
		psi_.set(v,someBasis);
	}

	const VectorWithOffsetType& gs() const { return psi_; }

	bool includeGroundStage() const {return true; }

	SizeType size() const
	{
		if (!correctionEnabled_) return 0;
		if (commonTargetting_.allStages(DISABLED,stage_)) return 1;
		return targetVectors_.size();
	}

	const VectorWithOffsetType& operator()(SizeType i) const
	{
		return targetVectors_[i];
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
		// //corner case
		SizeType x = (site==1) ? 0 : numberOfSites-1;
		evolve(Eg,direction,x,loopNumber);
	}

	void evolve(RealType Eg,SizeType direction,SizeType site,
	            SizeType loopNumber)
	{
		SizeType count =0;
		VectorWithOffsetType phiOld = psi_;
		VectorWithOffsetType phiNew;
		VectorWithOffsetType vectorSum;

		SizeType max = tstStruct_.sites.size();

		if (commonTargetting_.noStageIs(DISABLED,stage_)) max = 1;

		// Loop over each operator that needs to be applied
		// in turn to the g.s.
		for (SizeType i=0;i<max;i++) {
			count += evolve(i,phiNew,phiOld,Eg,direction,site,loopNumber,max-1);
			if (tstStruct_.concatenation==PRODUCT) {
				phiOld = phiNew;
			} else {
				if (i==0){
					vectorSum = phiNew;
				} else {
					if (tstStruct_.sites[i]==site)
						vectorSum += phiNew;
				}
			}
		}
		if (tstStruct_.concatenation==SUM) phiNew = vectorSum;

		if (!commonTargetting_.noStageIs(DISABLED,stage_))
			Eg_ = Eg;

		if (direction!=INFINITE) {
			correctionEnabled_=true;
			typename PsimagLite::Vector<SizeType>::Type block1(1,site);
			addCorrection(direction,block1);
		}

		if (count==0) return;

		calcDynVectors(phiNew,direction);
	}

	void initialGuess(VectorWithOffsetType& v,
	                  const typename PsimagLite::Vector<SizeType>::Type& block) const
	{
		commonTargetting_.initialGuess(v,wft_,block,psi_,stage_,weight_,targetVectors_);
	}

	const LeftRightSuperType& leftRightSuper() const { return lrs_; }

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          IoOutputType& io) const
	{
		if (block.size()!=1)
			throw PsimagLite::RuntimeError("CorrectionVectorTargetting only supports blocks of size 1\n");
		SizeType type = tstStruct_.type;
		int fermionSign = commonTargetting_.findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		typename PostProcType::ParametersType params;
		params.Eg = Eg_;
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		PostProcType cf(ab_,params);

		commonTargetting_.save(block,io,cf,targetVectors_);
		psi_.save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		for (SizeType i=0;i<stage_.size();i++) stage_[i] = CONVERGING;

		typename IoType::In io(f);

		commonTargetting_.load(io,targetVectors_);

		psi_.load(io,"PSI");
	}

	RealType time() const { return 0; }

	void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
	{}

	bool end() const { return false; }

private:

	SizeType evolve(SizeType i,
	                VectorWithOffsetType& phiNew,
	                VectorWithOffsetType& phiOld,
	                RealType Eg,
	                SizeType direction,
	                SizeType site,
	                SizeType loopNumber,
	                SizeType lastI)
	{
		if (tstStruct_.startingLoops[i]>loopNumber || direction==INFINITE) return 0;

		if (site != tstStruct_.sites[i] && stage_[i]==DISABLED) return 0;

		if (site == tstStruct_.sites[i] && stage_[i]==DISABLED) stage_[i]=OPERATOR;
		else stage_[i]=CONVERGING;
		if (stage_[i] == OPERATOR) commonTargetting_.checkOrder(i,stage_);

		PsimagLite::OstringStream msg;
		msg<<"Evolving, stage="<<commonTargetting_.getStage(i,stage_);
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
		SizeType numberOfSites = lrs_.super().block().size();
		if (stage_[i]==OPERATOR) {

			BorderEnumType corner = (tstStruct_.sites[i]==0 ||
			               tstStruct_.sites[i]==numberOfSites -1) ?
			            ApplyOperatorType::BORDER_YES : ApplyOperatorType::BORDER_NO;

			PsimagLite::OstringStream msg;
			msg<<"I'm applying a local operator now";
			progress_.printline(msg,std::cout);
			typename PsimagLite::Vector<SizeType>::Type electrons;
			commonTargetting_.findElectronsOfOneSite(electrons,site);
			FermionSign fs(lrs_.left(),electrons);
			applyOpLocal_(phiNew,phiOld,tstStruct_.aOperators[i],
			              fs,systemOrEnviron,corner);
		} else if (stage_[i]== CONVERGING) {
			if (site==0 || site==numberOfSites -1)  {
				// don't wft since we did it before
				phiNew = targetVectors_[1];
				return;
			}
			PsimagLite::OstringStream msg;
			msg<<"I'm calling the WFT now";
			progress_.printline(msg,std::cout);

			phiNew.populateSectors(lrs_.super());

			// OK, now that we got the partition number right, let's wft:
			typename PsimagLite::Vector<SizeType>::Type nk(1,model_.hilbertSize(site));
			wft_.setInitialVector(phiNew,targetVectors_[1],lrs_,nk);
			phiNew.collapseSectors();

		} else {
			throw PsimagLite::RuntimeError("computePhi\n");
		}
	}

	void calcDynVectors(const VectorWithOffsetType& phi,
	                    SizeType systemOrEnviron)
	{
		for (SizeType i=1;i<targetVectors_.size();i++)
			targetVectors_[i] = phi;

		for (SizeType i=0;i<phi.sectors();i++) {
			VectorType sv;
			SizeType i0 = phi.sector(i);
			phi.extract(sv,i0);
			// g.s. is included separately
			// set Aq
			targetVectors_[1].setDataInSector(sv,i0);
			// set xi
			SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
			VectorType xi(sv.size(),0),xr(sv.size(),0);
			computeXiAndXr(xi,xr,sv,p);
			targetVectors_[2].setDataInSector(xi,i0);
			//set xr
			targetVectors_[3].setDataInSector(xr,i0);
			DenseMatrixType V;
			getLanczosVectors(V,sv,p);
		}
		setWeights();
		weightForContinuedFraction_ = phi*phi;
	}

	void getLanczosVectors(DenseMatrixType& V,
	                       const VectorType& sv,
	                       SizeType p)
	{
		typename ModelType::ModelHelperType modelHelper(p,lrs_);
		typedef typename LanczosSolverType::LanczosMatrixType
		        LanczosMatrixType;
		LanczosMatrixType h(&model_,&modelHelper);

		ParametersForSolverType params;
		params.steps = tstStruct_.steps;
		params.tolerance = tstStruct_.eps;
		params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

		LanczosSolverType lanczosSolver(h,params,&V);

		lanczosSolver.decomposition(sv,ab_);
		//calcIntensity(Eg,sv,V,ab);
	}

	void computeXiAndXr(VectorType& xi,
	                    VectorType& xr,
	                    const VectorType& sv,
	                    SizeType p)
	{
		typename ModelType::ModelHelperType modelHelper(p,lrs_);
		LanczosMatrixType h(&model_,&modelHelper);
		CorrectionVectorFunctionType cvft(h,tstStruct_);

		cvft.getXi(xi,sv);
		// make sure xr is zero
		for (SizeType i=0;i<xr.size();i++) xr[i] = 0;
		h.matrixVectorProduct(xr,xi);
		xr -= tstStruct_.omega*xi;
		xr /= tstStruct_.eta;
	}

	void guessPhiSectors(VectorWithOffsetType& phi,SizeType i,SizeType systemOrEnviron)
	{
		FermionSign fs(lrs_.left(),tstStruct_.electrons);
		if (allStages(CONVERGING)) {
			VectorWithOffsetType tmpVector = psi_;
			for (SizeType j=0;j<tstStruct_.aOperators.size();j++) {
				applyOpLocal_(phi,tmpVector,tstStruct_.aOperators[j],fs,
				              systemOrEnviron);
				tmpVector = phi;
			}
			return;
		}
		applyOpLocal_(phi,psi_,tstStruct_.aOperators[i],fs,
		              systemOrEnviron);
	}

	void setWeights()
	{
		gsWeight_ = commonTargetting_.setGsWeight(0.5);

		RealType sum  = 0;
		weight_.resize(targetVectors_.size());
		for (SizeType r=1;r<weight_.size();r++) {
			weight_[r] =0;
			for (SizeType i=0;i<targetVectors_[1].sectors();i++) {
				VectorType v,w;
				SizeType i0 = targetVectors_[1].sector(i);
				targetVectors_[1].extract(v,i0);
				targetVectors_[r].extract(w,i0);
				weight_[r] += dynWeightOf(v,w);
			}
			sum += weight_[r];
		}
		for (SizeType r=0;r<weight_.size();r++) weight_[r] *= (1.0 - gsWeight_)/sum;
	}

	RealType dynWeightOf(VectorType& v,const VectorType& w) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<v.size();i++) {
			RealType tmp = std::real(v[i]*w[i]);
			sum += tmp*tmp;
		}
		return sum;
	}

	void addCorrection(SizeType direction,const BlockType& block1)
	{
		commonTargetting_.computeCorrection(targetVectors_[0],direction,block1,psi_);
		weight_.resize(1);
		weight_[0]=tstStruct_.correctionA;
		gsWeight_ = 1.0-weight_[0];
	}

	typename PsimagLite::Vector<SizeType>::Type stage_;
	VectorWithOffsetType psi_;
	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	ApplyOperatorType applyOpLocal_;
	RealType gsWeight_;
	typename PsimagLite::Vector<VectorWithOffsetType>::Type targetVectors_;
	CommonTargettingType commonTargetting_;
	bool correctionEnabled_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	TridiagonalMatrixType ab_;
	RealType Eg_;
	RealType weightForContinuedFraction_;

}; // class CorrectionVectorTargetting

template<template<typename,typename,typename> class LanczosSolverTemplate,
         template<typename,typename> class InternalProductTemplate,
         template<typename,typename> class WaveFunctionTransfTemplate,
         typename ModelType_,
         typename IoType_,
         template<typename> class VectorWithOffsetTemplate>
std::ostream& operator<<(std::ostream& os,
                         const CorrectionVectorTargetting<LanczosSolverTemplate,
                         InternalProductTemplate,
                         WaveFunctionTransfTemplate,ModelType_,IoType_,
                         VectorWithOffsetTemplate>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // CORRECTION_VECTOR_TARG_H

