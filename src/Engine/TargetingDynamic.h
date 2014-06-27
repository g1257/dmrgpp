
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

/*! \file TargetingDynamic.h
 *
 * Implements the targetting required by
 * a simple continued fraction calculation
 * of dynamical observables
 *
 */

#ifndef TARGETING_DYNAMIC_H
#define TARGETING_DYNAMIC_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "ParametersForSolver.h"
#include "TargetParamsDynamic.h"
#include "VectorWithOffsets.h"
#include "TargetingBase.h"
#include <cassert>
#include "Concurrency.h"
#include "Parallelizer.h"
#include "ProgramGlobals.h"
#include "ParallelWft.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_>
class TargetingDynamic : public TargetingBase<LanczosSolverTemplate,
                                              MatrixVectorType_,
                                              WaveFunctionTransfType_> {

public:

	typedef TargetingBase<LanczosSolverTemplate,
	                      MatrixVectorType_,
	                      WaveFunctionTransfType_> BaseType;
	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef TargetParamsDynamic<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              MatrixVectorType,
	                              VectorType> LanczosSolverType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;

	enum {DISABLED,OPERATOR,CONVERGING};
	enum {
		EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
		INFINITE=WaveFunctionTransfType::INFINITE
	};

	static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
	static SizeType const SUM = TargettingParamsType::SUM;

	TargetingDynamic(const LeftRightSuperType& lrs,
	                 const ModelType& model,
	                 const TargettingParamsType& tstStruct,
	                 const WaveFunctionTransfType& wft,
	                 const SizeType& quantumSector,
	                 InputValidatorType& io)
	    : BaseType(lrs,model,tstStruct,wft,0,0),
	      tstStruct_(tstStruct),
	      wft_(wft),
	      progress_("TargetingDynamic"),
	      gsWeight_(1.0),
	      paramsForSolver_(io,"DynamicDmrg"),
	      weightForContinuedFraction_(0)
	{
		if (!wft.isEnabled())
			throw PsimagLite::RuntimeError(" TargetingDynamic needs an enabled wft\n");
	}

	RealType weight(SizeType i) const
	{
		assert(!this->common().allStages(DISABLED));
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (this->common().allStages(DISABLED)) return 1.0;
		return gsWeight_;
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
		SizeType numberOfSites = this->leftRightSuper().super().block().size();
		if (site>1 && site<numberOfSites-2) return;
		// //corner case
		SizeType x = (site==1) ? 0 : numberOfSites-1;
		evolve(Eg,direction,x,loopNumber);
	}

	template<typename IoOutputType>
	void save(const typename PsimagLite::Vector<SizeType>::Type& block,
	          IoOutputType& io) const
	{
		assert(block.size()==1);

		SizeType type = tstStruct_.type();
		int fermionSign = this->common().findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		if (ab_.size()<2) return;
		typename PostProcType::ParametersType params = paramsForSolver_;
		params.Eg = this->common().energy();
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		if (tstStruct_.aOperators()[0].fermionSign>0) s2 *= s;

		PostProcType cf(ab_,reortho_,params);

		PsimagLite::String str = "#TCENTRALSITE=" + ttos(block[0]);
		io.printline(str);

		this->common().save(block,io,cf,this->common().targetVectors());

		this->common().psi().save(io,"PSI");
	}

	void load(const PsimagLite::String& f)
	{
		this->common().template load<TimeSerializerType>(f);
	}

private:

	void evolve(RealType Eg,
	            SizeType direction,
	            SizeType site,
	            SizeType loopNumber)
	{

		VectorWithOffsetType phiNew;
		SizeType count = this->common().getPhi(phiNew,Eg,direction,site,loopNumber);

		if (count==0) return;

		calcLanczosVectors(gsWeight_,weight_,phiNew,direction);

		typename PsimagLite::Vector<SizeType>::Type block(1,site);
		this->common().cocoon(block,direction);
	}

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
			SizeType p = this->leftRightSuper().super().findPartitionNumber(phi.offset(i0));
			getLanczosVectors(V,sv,p);
			if (i==0) {
				assert(V.n_col() > 0);
				this->common().targetVectorsResize(V.n_col());
				for (SizeType j=0;j<this->common().targetVectors().size();j++)
					this->common().targetVectors(j) = phi;
			}
			setVectors(V,i0);
		}

		setWeights();
		if (fabs(weightForContinuedFraction_)<1e-6)
			weightForContinuedFraction_ = std::real(phi*phi);
	}

	void wftLanczosVectors(SizeType site,const VectorWithOffsetType& phi)
	{
		this->common().targetVectors()[0] = phi;
		// don't wft since we did it before
		SizeType numberOfSites = this->leftRightSuper().super().block().size();
		if (site==0 || site==numberOfSites -1)  return;

		typedef ParallelWft<VectorWithOffsetType,
		                    WaveFunctionTransfType,
		                    LeftRightSuperType> ParallelWftType;
		typedef PsimagLite::Parallelizer<ParallelWftType> ParallelizerType;
		ParallelizerType threadedWft(PsimagLite::Concurrency::npthreads,
		                             PsimagLite::MPI::COMM_WORLD);

		ParallelWftType helperWft(this->common().targetVectors(),this->model().hilbertSize(site),wft_,this->leftRightSuper());
		threadedWft.loopCreate(this->common().targetVectors().size()-1,helperWft,this->model().concurrency());

		for (SizeType i=1;i<this->common().targetVectors().size();i++) {
			assert(this->common().targetVectors()[i].size()==this->common().targetVectors()[0].size());
		}
	}

	void getLanczosVectors(DenseMatrixType& V,
	                       const VectorType& sv,
	                       SizeType p)
	{
		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,this->leftRightSuper(),threadId);
		typename LanczosSolverType::LanczosMatrixType h(&this->model(),&modelHelper);
		paramsForSolver_.lotaMemory = true;
		LanczosSolverType lanczosSolver(h,paramsForSolver_,&V);

		lanczosSolver.decomposition(sv,ab_);

		reortho_ = lanczosSolver.reorthogonalizationMatrix();
	}

	void setVectors(const DenseMatrixType& V,
	                SizeType i0)
	{
		for (SizeType i=0;i<this->common().targetVectors().size();i++) {
			VectorType tmp(V.n_row());
			for (SizeType j=0;j<tmp.size();j++) tmp[j] = V(j,i);
			this->common().targetVectors(i).setDataInSector(tmp,i0);
		}
	}

	void setWeights()
	{
		gsWeight_ = 0.0;
		RealType sum  = 0;
		weight_.resize(this->common().targetVectors().size());
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

	const TargettingParamsType& tstStruct_;
	const WaveFunctionTransfType& wft_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	ParametersForSolverType paramsForSolver_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	TridiagonalMatrixType ab_;
	DenseMatrixRealType reortho_;
	RealType weightForContinuedFraction_;
}; // class TargetingDynamic

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingDynamic<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // TARGETING_DYNAMIC_H

