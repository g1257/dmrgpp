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
#include "CorrectionVectorParams.h"
#include "VectorWithOffsets.h"
#include "CorrectionVectorFunction.h"
#include "CommonTargetting.h"
#include "ParametersForSolver.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_,
         typename IoType_>
class CorrectionVectorTargetting  {

public:

	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef IoType_ IoType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType
	LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef CorrectionVectorParams<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              MatrixVectorType,
	                              VectorType> LanczosSolverType;
	typedef VectorType TargetVectorType;
	typedef typename IoType::In IoInputType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef DynamicSerializer<VectorWithOffsetType,PostProcType> DynamicSerializerType;
	typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,
	                                 TargettingParamsType> CorrectionVectorFunctionType;
	typedef TargetHelper<ModelType,
	                     TargettingParamsType,
	                     WaveFunctionTransfType> TargetHelperType;
	typedef CommonTargetting<TargetHelperType,
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
	    : lrs_(lrs),
	      model_(model),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      progress_("CorrectionVectorTargetting"),
	      gsWeight_(1.0),
	      commonTargetting_(lrs,model,tstStruct,wft,4),
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
		return commonTargetting_.normSquared(i);
	}

	template<typename SomeBasisType>
	void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
	           const SomeBasisType& someBasis)
	{
		commonTargetting_.psi().set(v,someBasis);
	}

	const VectorWithOffsetType& gs() const { return commonTargetting_.psi(); }

	bool includeGroundStage() const {return true; }

	SizeType size() const
	{
		if (!correctionEnabled_) return 0;
		if (commonTargetting_.allStages(DISABLED)) return 1;
		return commonTargetting_.targetVectors().size();
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
		// //corner case
		SizeType x = (site==1) ? 0 : numberOfSites-1;
		evolve(Eg,direction,x,loopNumber);
	}

	void evolve(RealType Eg,SizeType direction,SizeType site,
	            SizeType loopNumber)
	{
		VectorWithOffsetType phiNew;
		SizeType count = commonTargetting_.getPhi(phiNew,Eg,direction,site,loopNumber);

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
		commonTargetting_.initialGuess(v,block);
	}

	const LeftRightSuperType& leftRightSuper() const { return lrs_; }

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          IoOutputType& io) const
	{
		if (block.size()!=1)
			throw PsimagLite::RuntimeError("CorrectionVectorTargetting only supports blocks of size 1\n");
		SizeType type = tstStruct_.type();
		int fermionSign = commonTargetting_.findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		typename PostProcType::ParametersType params;
		params.Eg = commonTargetting_.energy();
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		PostProcType cf(ab_,reortho_,params);

		commonTargetting_.save(block,io,cf,commonTargetting_.targetVectors());
		commonTargetting_.psi().save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		commonTargetting_.setAllStagesTo(CONVERGING);

		IoInputType io(f);

		TimeSerializerType dynS(io,IoInputType::LAST_INSTANCE);
		commonTargetting_.loadTargetVectors(dynS);

		commonTargetting_.psi().load(io,"PSI");
	}

	RealType time() const { return 0; }

	void updateOnSiteForTimeDep(BasisWithOperatorsType& basisWithOps) const
	{}

	bool end() const { return false; }

private:

	void calcDynVectors(const VectorWithOffsetType& phi,
	                    SizeType systemOrEnviron)
	{
		for (SizeType i=1;i<commonTargetting_.targetVectors().size();i++)
			commonTargetting_.targetVectors(i) = phi;

		for (SizeType i=0;i<phi.sectors();i++) {
			VectorType sv;
			SizeType i0 = phi.sector(i);
			phi.extract(sv,i0);
			// g.s. is included separately
			// set Aq
			commonTargetting_.targetVectors(1).setDataInSector(sv,i0);
			// set xi
			SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
			VectorType xi(sv.size(),0),xr(sv.size(),0);
			computeXiAndXr(xi,xr,sv,p);
			commonTargetting_.targetVectors(2).setDataInSector(xi,i0);
			//set xr
			commonTargetting_.targetVectors(3).setDataInSector(xr,i0);
			DenseMatrixType V;
			getLanczosVectors(V,sv,p);
		}
		setWeights();
		weightForContinuedFraction_ = std::real(phi*phi);
	}

	void getLanczosVectors(DenseMatrixType& V,
	                       const VectorType& sv,
	                       SizeType p)
	{
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,threadId);
		typedef typename LanczosSolverType::LanczosMatrixType
		        LanczosMatrixType;
		LanczosMatrixType h(&model_,&modelHelper);

		ParametersForSolverType params;
		params.steps = tstStruct_.steps();
		params.tolerance = tstStruct_.eps();
		params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

		LanczosSolverType lanczosSolver(h,params,&V);

		lanczosSolver.decomposition(sv,ab_);
		reortho_ = lanczosSolver.reorthogonalizationMatrix();
	}

	void computeXiAndXr(VectorType& xi,
	                    VectorType& xr,
	                    const VectorType& sv,
	                    SizeType p)
	{
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,threadId);
		LanczosMatrixType h(&model_,&modelHelper);
		CorrectionVectorFunctionType cvft(h,tstStruct_);

		cvft.getXi(xi,sv);
		// make sure xr is zero
		for (SizeType i=0;i<xr.size();i++) xr[i] = 0;
		h.matrixVectorProduct(xr,xi);
		xr -= tstStruct_.omega()*xi;
		xr /= tstStruct_.eta();
	}

	void setWeights()
	{
		gsWeight_ = commonTargetting_.setGsWeight(0.5);

		RealType sum  = 0;
		weight_.resize(commonTargetting_.targetVectors().size());
		for (SizeType r=1;r<weight_.size();r++) {
			weight_[r] =0;
			for (SizeType i=0;i<commonTargetting_.targetVectors()[1].sectors();i++) {
				VectorType v,w;
				SizeType i0 = commonTargetting_.targetVectors()[1].sector(i);
				commonTargetting_.targetVectors()[1].extract(v,i0);
				commonTargetting_.targetVectors()[r].extract(w,i0);
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
		commonTargetting_.computeCorrection(direction,block1);
		weight_.resize(1);
		weight_[0]=tstStruct_.correctionA();
		gsWeight_ = 1.0-weight_[0];
	}

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	CommonTargettingType commonTargetting_;
	bool correctionEnabled_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	TridiagonalMatrixType ab_;
	DenseMatrixRealType reortho_;
	RealType weightForContinuedFraction_;

}; // class CorrectionVectorTargetting

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType,
         typename IoType_>
std::ostream& operator<<(std::ostream& os,
                         const CorrectionVectorTargetting<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType,IoType_>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // CORRECTION_VECTOR_TARG_H

