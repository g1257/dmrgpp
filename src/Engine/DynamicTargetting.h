
/*
Copyright (c) 2009, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file DynamicTargetting.h
 *
 * Implements the targetting required by
 * a simple continued fraction calculation
 * of dynamical observables
 *
 */

#ifndef DYNAMICTARGETTING_H
#define DYNAMICTARGETTING_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ParametersForSolver.h"
#include "DynamicDmrgParams.h"
#include "VectorWithOffsets.h"
#include "CommonTargetting.h"
#include <cassert>
#include "Concurrency.h"
#include "Parallelizer.h"
#include "ProgramGlobals.h"
#include "ParallelWft.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_,
         typename IoType_>
class DynamicTargetting  {

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
	typedef DynamicDmrgParams<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef typename IoType::In IoInputType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              MatrixVectorType,
	                              VectorType> LanczosSolverType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef TargetHelper<ModelType,
	                     TargettingParamsType,
	                     WaveFunctionTransfType> TargetHelperType;
	typedef CommonTargetting<TargetHelperType,
	                         VectorWithOffsetType,
	                         LanczosSolverType> CommonTargettingType;

	enum {DISABLED,OPERATOR,CONVERGING};
	enum {
		EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
		INFINITE=WaveFunctionTransfType::INFINITE
	};

	static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
	static SizeType const SUM = TargettingParamsType::SUM;

	DynamicTargetting(const LeftRightSuperType& lrs,
	                  const ModelType& model,
	                  const TargettingParamsType& tstStruct,
	                  const WaveFunctionTransfType& wft,
	                  const SizeType& quantumSector)
	    : lrs_(lrs),
	      model_(model),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      progress_("DynamicTargetting"),
	      gsWeight_(1.0),
	      commonTargetting_(lrs,model,tstStruct,wft),
	      weightForContinuedFraction_(0)
	{
		if (!wft.isEnabled())
			throw PsimagLite::RuntimeError(" DynamicTargetting needs an enabled wft\n");

		paramsForSolver_.steps = tstStruct_.steps();
		paramsForSolver_.tolerance = tstStruct_.eps();
		paramsForSolver_.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;
	}

	RealType weight(SizeType i) const
	{
		assert(!commonTargetting_.allStages(DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (commonTargetting_.allStages(DISABLED)) return 1.0;
		return gsWeight_;
	}

	template<typename SomeBasisType>
	void setGs(const typename PsimagLite::Vector<TargetVectorType>::Type& v,
	           const SomeBasisType& someBasis)
	{
		commonTargetting_.psi().set(v,someBasis);
	}

	RealType normSquared(SizeType i) const
	{
		return commonTargetting_.normSquared(i);
	}

	const VectorWithOffsetType& gs() const
	{
		return commonTargetting_.psi();
	}

	bool includeGroundStage() const {return fabs(gsWeight_)>1e-6; }

	SizeType size() const
	{
		if (commonTargetting_.allStages(DISABLED)) return 0;
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

	void evolve(RealType Eg,
	            SizeType direction,
	            SizeType site,
	            SizeType loopNumber)
	{

		VectorWithOffsetType phiNew;
		SizeType count = commonTargetting_.getPhi(phiNew,Eg,direction,site,loopNumber);

		if (count==0) return;

		calcLanczosVectors(gsWeight_,weight_,phiNew,direction);

		if (model_.params().insitu=="" || !includeGroundStage()) return;

		if (BasisType::useSu2Symmetry()) {
			commonTargetting_.noCocoon("not when SU(2) symmetry is in use");
			return;
		}

		try {
			commonTargetting_.cocoon(direction,site,commonTargetting_.psi(),"PSI");
		} catch (std::exception& e) {
			commonTargetting_.noCocoon("unsupported by the model");
		}
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
		assert(block.size()==1);

		SizeType type = tstStruct_.type();
		int fermionSign = commonTargetting_.findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		if (ab_.size()<2) return;
		typename PostProcType::ParametersType params = paramsForSolver_;
		params.Eg = commonTargetting_.energy();
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		if (tstStruct_.aOperators()[0].fermionSign>0) s2 *= s;

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
			commonTargetting_.setAllStagesTo(CONVERGING);
			TimeSerializerType dynS(io,IoInputType::LAST_INSTANCE);
			commonTargetting_.loadTargetVectors(dynS);
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

	const ModelType& model() const { return model_; }

	bool end() const { return false; }

private:

	void calcLanczosVectors(RealType& gsWeight,
	                        typename PsimagLite::Vector<RealType>::Type& weights,
	                        const VectorWithOffsetType& phi,
	                        SizeType systemOrEnviron)
	{
		for (SizeType i=0;i<phi.sectors();i++) {
			VectorType sv;
			SizeType i0 = phi.sector(i);
			phi.extract(sv,i0);
			DenseMatrixType V;
			SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
			getLanczosVectors(V,sv,p);
			if (i==0) {
				commonTargetting_.targetVectorsResize(V.n_col());
				for (SizeType j=0;j<commonTargetting_.targetVectors().size();j++)
					commonTargetting_.targetVectors(j) = phi;
			}
			setVectors(V,i0);
		}

		setWeights();
		if (fabs(weightForContinuedFraction_)<1e-6)
			weightForContinuedFraction_ = std::real(phi*phi);
		//			std::cerr<<"weight==============="<<weightForContinuedFraction_<<"\n";
	}

	void wftLanczosVectors(SizeType site,const VectorWithOffsetType& phi)
	{
		commonTargetting_.targetVectors()[0] = phi;
		// don't wft since we did it before
		SizeType numberOfSites = lrs_.super().block().size();
		if (site==0 || site==numberOfSites -1)  return;

		typedef ParallelWft<VectorWithOffsetType,
		                    WaveFunctionTransfType,
		                    LeftRightSuperType> ParallelWftType;
		typedef PsimagLite::Parallelizer<ParallelWftType> ParallelizerType;
		ParallelizerType threadedWft(PsimagLite::Concurrency::npthreads,
		                             PsimagLite::MPI::COMM_WORLD);

		ParallelWftType helperWft(commonTargetting_.targetVectors(),model_.hilbertSize(site),wft_,lrs_);
		threadedWft.loopCreate(commonTargetting_.targetVectors().size()-1,helperWft,model_.concurrency());

		for (SizeType i=1;i<commonTargetting_.targetVectors().size();i++) {
			assert(commonTargetting_.targetVectors()[i].size()==commonTargetting_.targetVectors()[0].size());
		}
	}

	void getLanczosVectors(DenseMatrixType& V,
	                       const VectorType& sv,
	                       SizeType p)
	{
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,threadId);
		typename LanczosSolverType::LanczosMatrixType h(&model_,&modelHelper);

		LanczosSolverType lanczosSolver(h,paramsForSolver_,&V);

		lanczosSolver.decomposition(sv,ab_);

		reortho_ = lanczosSolver.reorthogonalizationMatrix();
	}

	void setVectors(const DenseMatrixType& V,
	                SizeType i0)
	{
		for (SizeType i=0;i<commonTargetting_.targetVectors().size();i++) {
			VectorType tmp(V.n_row());
			for (SizeType j=0;j<tmp.size();j++) tmp[j] = V(j,i);
			commonTargetting_.targetVectors(i).setDataInSector(tmp,i0);
		}
	}

	void setWeights()
	{
		gsWeight_ = 0.0;
		RealType sum  = 0;
		weight_.resize(commonTargetting_.targetVectors().size());
		for (SizeType r=0;r<weight_.size();r++) {
			weight_[r] = 1.0;
			sum += weight_[r];
		}

		for (SizeType r=0;r<weight_.size();r++) weight_[r] *=(1.0-gsWeight_)/sum;
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

	const LeftRightSuperType& lrs_;
	const ModelType& model_;
	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	CommonTargettingType commonTargetting_;
	ParametersForSolverType paramsForSolver_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	TridiagonalMatrixType ab_;
	DenseMatrixRealType reortho_;
	RealType weightForContinuedFraction_;
}; // class DynamicTargetting

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType,
         typename IoType_>
std::ostream& operator<<(std::ostream& os,
                         const DynamicTargetting<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType,IoType_>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // DYNAMICTARGETTING_H

