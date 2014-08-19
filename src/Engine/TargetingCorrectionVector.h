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

/*! \file TargetingCorrectionVector.h
 *
 * Implements the targetting required by
 * the correction targetting method
 *
 */

#ifndef TARGETING_CORRECTION_VECTOR_H
#define TARGETING_CORRECTION_VECTOR_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "TargetParamsCorrectionVector.h"
#include "VectorWithOffsets.h"
#include "CorrectionVectorFunction.h"
#include "TargetingBase.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "TimeSerializer.h"
#include "FreqEnum.h"

namespace Dmrg {

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType_,
         typename WaveFunctionTransfType_>
class TargetingCorrectionVector : public TargetingBase<LanczosSolverTemplate,
                                                       MatrixVectorType_,
                                                       WaveFunctionTransfType_> {

	typedef TargetingBase<LanczosSolverTemplate,
	                      MatrixVectorType_,
	                      WaveFunctionTransfType_> BaseType;

	class CalcR {

		typedef typename MatrixVectorType_::ModelType ModelType;
		typedef typename ModelType::RealType RealType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
		typedef TargetParamsCorrectionVector<ModelType> TargettingParamsType;

		class Action {

		public:

			enum ActionEnum {ACTION_IMAG, ACTION_REAL};

			Action(const TargettingParamsType& tstStruct,
			       RealType E0,
			       const VectorRealType& eigs)
			    : tstStruct_(tstStruct),E0_(E0),eigs_(eigs)
			{}

			RealType operator()(SizeType k) const
			{
				 if (tstStruct_.omega().first == PsimagLite::FREQ_REAL)
					 return actionWhenReal(k);

				 return actionWhenMatsubara(k);
			}

			void setReal() const
			{
				action_ = ACTION_REAL;
			}

			void setImag() const
			{
				action_ = ACTION_IMAG;
			}

		private:

			RealType actionWhenReal(SizeType k) const
			{
				RealType sign = (tstStruct_.type() == 0) ? -1.0 : 1.0;
				RealType part1 =  (eigs_[k] - E0_)*sign + tstStruct_.omega().second;
				RealType denom = part1*part1 + tstStruct_.eta()*tstStruct_.eta();
				return (action_ == ACTION_IMAG) ? tstStruct_.eta()/denom : -part1 / denom;
			}

			RealType actionWhenMatsubara(SizeType k) const
			{
				RealType sign = (tstStruct_.type() == 0) ? -1.0 : 1.0;
				RealType wn = tstStruct_.omega().second;
				RealType part1 =  (eigs_[k] - E0_)*sign;
				RealType denom = part1*part1 + wn*wn;
				return (action_ == ACTION_IMAG) ? wn/denom : -part1 / denom;
			}

			const TargettingParamsType& tstStruct_;
			RealType E0_;
			const VectorRealType& eigs_;
			mutable ActionEnum action_;
		};

	public:

		typedef Action ActionType;

		CalcR(const TargettingParamsType& tstStruct,
		      RealType E0,
		      const VectorRealType& eigs)
		    : action_(tstStruct,E0,eigs)
		{}

		const Action& imag() const
		{
			action_.setImag();
			return action_;
		}

		const Action& real() const
		{
			action_.setReal();
			return action_;
		}

	private:

		Action action_;
	};

	typedef CalcR CalcRType;

public:

	typedef MatrixVectorType_ MatrixVectorType;
	typedef typename MatrixVectorType::ModelType ModelType;
	typedef typename ModelType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename ModelType::OperatorsType OperatorsType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef TargetParamsCorrectionVector<ModelType> TargettingParamsType;
	typedef typename BasisType::BlockType BlockType;
	typedef WaveFunctionTransfType_ WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::ParametersForSolver<RealType> ParametersForSolverType;
	typedef LanczosSolverTemplate<ParametersForSolverType,
	                              MatrixVectorType,
	                              VectorType> LanczosSolverType;
	typedef VectorType TargetVectorType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef DynamicSerializer<VectorWithOffsetType,PostProcType> DynamicSerializerType;
	typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,
	                                 TargettingParamsType> CorrectionVectorFunctionType;
	typedef ParallelTriDiag<ModelType,
	                        LanczosSolverType,
	                        VectorWithOffsetType> ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;

	enum {DISABLED,OPERATOR,CONVERGING};

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
	      EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
	      INFINITE=WaveFunctionTransfType::INFINITE};

	static SizeType const PRODUCT = TargettingParamsType::PRODUCT;
	static SizeType const SUM = TargettingParamsType::SUM;

	TargetingCorrectionVector(const LeftRightSuperType& lrs,
	                          const ModelType& model,
	                          const TargettingParamsType& tstStruct,
	                          const WaveFunctionTransfType& wft,
	                          const SizeType& quantumSector,
	                          InputValidatorType& ioIn)
	    : BaseType(lrs,model,tstStruct,wft,4,1),
	      tstStruct_(tstStruct),
	      ioIn_(ioIn),
	      progress_("TargetingCorrectionVector"),
	      gsWeight_(1.0),
	      correctionEnabled_(false),
	      paramsForSolver_(ioIn,"DynamicDmrg")
	{
		if (!wft.isEnabled())
			throw PsimagLite::RuntimeError("TargetingCorrectionVector needs wft\n");
	}

	RealType weight(SizeType i) const
	{
		return weight_[i];
	}

	RealType gsWeight() const
	{
		if (!correctionEnabled_) return 1.0;
		return gsWeight_;
	}

	SizeType size() const
	{
		if (!correctionEnabled_) return 0;
		return BaseType::size();
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
		if (block.size()!=1)
			throw PsimagLite::RuntimeError("TargetingCorrectionVector only supports blocks of size 1\n");
		SizeType type = tstStruct_.type();
		int fermionSign = this->common().findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		typename PostProcType::ParametersType params(ioIn_,"DynamicDmrg");
		params.Eg = this->common().energy();
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		PostProcType cf(ab_,reortho_,params);

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

		if (direction!=INFINITE) {
			correctionEnabled_=true;
			typename PsimagLite::Vector<SizeType>::Type block1(1,site);
			addCorrection(direction,block1);
		}

		if (count==0) return;

		calcDynVectors(phiNew,direction);
		for (SizeType i = 1; i < this->common().targetVectors().size(); ++i) {
			PsimagLite::String label = "P" + ttos(i);
			this->common().cocoon(direction,
			                      site,
			                      this->common().psi(),
			                      "PSI",
			                      this->common().targetVectors(i),
			                      label);
		}
	}

	void calcDynVectors(const VectorWithOffsetType& phi,
	                    SizeType systemOrEnviron)
	{
		for (SizeType i=1;i<this->common().targetVectors().size();i++)
			this->common().targetVectors(i) = phi;

		VectorMatrixFieldType V(phi.sectors());
		VectorMatrixFieldType T(phi.sectors());

		VectorSizeType steps(phi.sectors());

		triDiag(phi,T,V,steps);

		VectorVectorRealType eigs(phi.sectors());

		for (SizeType ii=0;ii<phi.sectors();ii++)
			PsimagLite::diag(T[ii],eigs[ii],'V');

		for (SizeType i=0;i<phi.sectors();i++) {
			VectorType sv;
			SizeType i0 = phi.sector(i);
			phi.extract(sv,i0);
			// g.s. is included separately
			// set Aq
			this->common().targetVectors(1).setDataInSector(sv,i0);
			// set xi
			SizeType p = this->leftRightSuper().super().findPartitionNumber(phi.offset(i0));
			VectorType xi(sv.size(),0),xr(sv.size(),0);

			if (tstStruct_.algorithm() == TargettingParamsType::KRYLOV) {
				computeXiAndXrKrylov(xi,xr,phi,i0,V[i],T[i],eigs[i],steps[i]);
			} else {
				computeXiAndXrIndirect(xi,xr,sv,p);
			}

			this->common().targetVectors(2).setDataInSector(xi,i0);
			//set xr
			this->common().targetVectors(3).setDataInSector(xr,i0);
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
		typename ModelType::ModelHelperType modelHelper(p,this->leftRightSuper(),threadId);
		typedef typename LanczosSolverType::LanczosMatrixType
		        LanczosMatrixType;
		LanczosMatrixType h(&this->model(),&modelHelper);
		LanczosSolverType lanczosSolver(h,paramsForSolver_,&V);

		lanczosSolver.decomposition(sv,ab_);
		reortho_ = lanczosSolver.reorthogonalizationMatrix();
	}

	void computeXiAndXrIndirect(VectorType& xi,
	                            VectorType& xr,
	                            const VectorType& sv,
	                            SizeType p)
	{
		if (tstStruct_.omega().first != PsimagLite::FREQ_REAL)
			throw PsimagLite::RuntimeError("Matsubara only with KRYLOV\n");

		SizeType threadId = 0;
		typename ModelType::ModelHelperType modelHelper(p,this->leftRightSuper(),threadId);
		LanczosMatrixType h(&this->model(),&modelHelper);
		RealType E0 = this->common().energy();
		CorrectionVectorFunctionType cvft(h,tstStruct_,E0);

		cvft.getXi(xi,sv);
		// make sure xr is zero
		for (SizeType i=0;i<xr.size();i++) xr[i] = 0;
		h.matrixVectorProduct(xr,xi);
		xr -= (tstStruct_.omega().second+E0)*xi;
		xr /= tstStruct_.eta();
	}

	void computeXiAndXrKrylov(VectorType& xi,
	                          VectorType& xr,
	                          const VectorWithOffsetType& phi,
	                          SizeType i0,
	                          const MatrixComplexOrRealType& V,
	                          const MatrixComplexOrRealType& T,
	                          const VectorRealType& eigs,
	                          SizeType steps)
	{
		SizeType n2 = steps;
		SizeType n = V.n_row();
		if (T.n_col()!=T.n_row()) throw PsimagLite::RuntimeError("T is not square\n");
		if (V.n_col()!=T.n_col()) throw PsimagLite::RuntimeError("V is not nxn2\n");
		// for (SizeType j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
		ComplexOrRealType zone = 1.0;
		ComplexOrRealType zzero = 0.0;

		TargetVectorType tmp(n2);
		VectorType r(n2);
		CalcRType what(tstStruct_,this->common().energy(),eigs);

		calcR(r,what.imag(),T,V,phi,eigs,n2,i0);

		psimag::BLAS::GEMV('N', n2, n2, zone, &(T(0,0)), n2, &(r[0]), 1, zzero, &(tmp[0]), 1 );

		xi.resize(n);
		psimag::BLAS::GEMV('N', n,  n2, zone, &(V(0,0)), n, &(tmp[0]),1, zzero, &(xi[0]),   1 );

		calcR(r,what.real(),T,V,phi,eigs,n2,i0);

		psimag::BLAS::GEMV('N', n2, n2, zone, &(T(0,0)), n2, &(r[0]), 1, zzero, &(tmp[0]), 1 );

		xr.resize(n);
		psimag::BLAS::GEMV('N', n,  n2, zone, &(V(0,0)), n, &(tmp[0]),1, zzero, &(xr[0]),   1 );
	}

	void calcR(TargetVectorType& r,
	           const typename CalcRType::ActionType& whatRorI,
			   const MatrixComplexOrRealType& T,
			   const MatrixComplexOrRealType& V,
			   const VectorWithOffsetType& phi,
			   const VectorRealType& eigs,
			   SizeType n2,
			   SizeType i0)
	{
		for (SizeType k=0;k<n2;k++) {
			ComplexOrRealType sum = 0.0;
			for (SizeType kprime=0;kprime<n2;kprime++) {
				ComplexOrRealType tmpV = calcVTimesPhi(kprime,V,phi,i0);
				sum += std::conj(T(kprime,k))*tmpV;
			}

			r[k] = sum * whatRorI(k);
		}
	}

	ComplexOrRealType calcVTimesPhi(SizeType kprime,
	                                const MatrixComplexOrRealType& V,
	                                const VectorWithOffsetType& phi,
	                                SizeType i0) const
	{
		ComplexOrRealType ret = 0;
		SizeType total = phi.effectiveSize(i0);

		for (SizeType j=0;j<total;j++)
			ret += std::conj(V(j,kprime))*phi.fastAccess(i0,j);
		return ret;
	}

	void triDiag(const VectorWithOffsetType& phi,
	             VectorMatrixFieldType& T,
	             VectorMatrixFieldType& V,
	             typename PsimagLite::Vector<SizeType>::Type& steps)
	{
		PsimagLite::String options = this->model().params().options;
		bool cTridiag = (options.find("concurrenttridiag") != PsimagLite::String::npos);

		if (cTridiag) {
			typedef PsimagLite::Parallelizer<ParallelTriDiagType> ParallelizerType;
			ParallelizerType threadedTriDiag(PsimagLite::Concurrency::npthreads,
			                                 PsimagLite::MPI::COMM_WORLD);

			ParallelTriDiagType helperTriDiag(phi,T,V,steps,this->leftRightSuper(),this->model(),ioIn_);

			threadedTriDiag.loopCreate(phi.sectors(),helperTriDiag);
		} else {
			typedef PsimagLite::NoPthreads<ParallelTriDiagType> ParallelizerType;
			ParallelizerType threadedTriDiag(1,0);

			ParallelTriDiagType helperTriDiag(phi,T,V,steps,this->leftRightSuper(),this->model(),ioIn_);

			threadedTriDiag.loopCreate(phi.sectors(),helperTriDiag);
		}
	}

	void setWeights()
	{
		gsWeight_ = this->common().setGsWeight(0.5);

		RealType sum  = 0;
		weight_.resize(this->common().targetVectors().size());
		for (SizeType r=1;r<weight_.size();r++) {
			weight_[r] =0;
			for (SizeType i=0;i<this->common().targetVectors()[1].sectors();i++) {
				VectorType v,w;
				SizeType i0 = this->common().targetVectors()[1].sector(i);
				this->common().targetVectors()[1].extract(v,i0);
				this->common().targetVectors()[r].extract(w,i0);
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
		weight_.resize(1);
		weight_[0]=tstStruct_.correctionA();
		this->common().computeCorrection(direction,block1);
		gsWeight_ = 1.0-weight_[0];
	}

	const TargettingParamsType& tstStruct_;
	InputValidatorType& ioIn_;
	PsimagLite::ProgressIndicator progress_;
	RealType gsWeight_;
	bool correctionEnabled_;
	typename PsimagLite::Vector<RealType>::Type weight_;
	TridiagonalMatrixType ab_;
	DenseMatrixRealType reortho_;
	RealType weightForContinuedFraction_;
	ParametersForSolverType paramsForSolver_;
}; // class TargetingCorrectionVector

template<template<typename,typename,typename> class LanczosSolverTemplate,
         typename MatrixVectorType,
         typename WaveFunctionTransfType>
std::ostream& operator<<(std::ostream& os,
                         const TargetingCorrectionVector<LanczosSolverTemplate,
                         MatrixVectorType,
                         WaveFunctionTransfType>& tst)
{
	os<<"DT=NothingToSeeHereYet\n";
	return os;
}

} // namespace
/*@}*/
#endif // TARGETING_CORRECTION_VECTOR_H

