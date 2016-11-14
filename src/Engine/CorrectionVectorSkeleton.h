/*
Copyright (c) 2009-2016, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file CorrectionVectorSkeleton.h
 *
 * Common code to TargetingCorrectionVector TargetingRixs[Static and Dynamic]
 *
 */

#ifndef CORRECTION_VECTOR_SKEL_H
#define CORRECTION_VECTOR_SKEL_H

#include "ProgressIndicator.h"
#include "BLAS.h"
#include "VectorWithOffsets.h"
#include "CorrectionVectorFunction.h"
#include "ParametersForSolver.h"
#include "ParallelTriDiag.h"
#include "TimeSerializer.h"
#include "FreqEnum.h"

namespace Dmrg {

template<typename LanczosSolverType_,
         typename VectorWithOffsetType_,
         typename TargetingBaseType,
         typename TargetParamsType>
class CorrectionVectorSkeleton {

	typedef LanczosSolverType_ LanczosSolverType;

	class CalcR {

		typedef typename LanczosSolverType::LanczosMatrixType::ModelType ModelType;
		typedef typename ModelType::RealType RealType;
		typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

		class Action {

		public:

			enum ActionEnum {ACTION_IMAG, ACTION_REAL};

			Action(const TargetParamsType& tstStruct,
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
				return (action_ == ACTION_IMAG) ? tstStruct_.eta()/denom :
				                                  -part1/denom;
			}

			RealType actionWhenMatsubara(SizeType k) const
			{
				RealType sign = (tstStruct_.type() == 0) ? -1.0 : 1.0;
				RealType wn = tstStruct_.omega().second;
				RealType part1 =  (eigs_[k] - E0_)*sign;
				RealType denom = part1*part1 + wn*wn;
				return (action_ == ACTION_IMAG) ? wn/denom : -part1 / denom;
			}

			const TargetParamsType& tstStruct_;
			RealType E0_;
			const VectorRealType& eigs_;
			mutable ActionEnum action_;
		};

	public:

		typedef Action ActionType;

		CalcR(const TargetParamsType& tstStruct,
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

	typedef typename TargetingBaseType::MatrixVectorType MatrixVectorType;
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
	typedef typename BasisType::BlockType BlockType;
	typedef typename TargetingBaseType::WaveFunctionTransfType WaveFunctionTransfType;
	typedef typename WaveFunctionTransfType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef VectorType TargetVectorType;
	typedef TimeSerializer<VectorWithOffsetType> TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef DynamicSerializer<VectorWithOffsetType,PostProcType> DynamicSerializerType;
	typedef typename LanczosSolverType::LanczosMatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,
	TargetParamsType> CorrectionVectorFunctionType;
	typedef ParallelTriDiag<ModelType,
	LanczosSolverType,
	VectorWithOffsetType> ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename TargetingBaseType::InputSimpleOutType InputSimpleOutType;

	enum {EXPAND_ENVIRON=WaveFunctionTransfType::EXPAND_ENVIRON,
		  EXPAND_SYSTEM=WaveFunctionTransfType::EXPAND_SYSTEM,
		  INFINITE=WaveFunctionTransfType::INFINITE};

	static SizeType const PRODUCT = TargetParamsType::PRODUCT;
	static SizeType const SUM = TargetParamsType::SUM;

	CorrectionVectorSkeleton(InputValidatorType& ioIn,
	                         const TargetParamsType& tstStruct,
	                         const ModelType& model,
	                         const LeftRightSuperType& lrs,
	                         RealType energy)
	    : ioIn_(ioIn),
	      tstStruct_(tstStruct),
	      model_(model),
	      lrs_(lrs),
	      energy_(energy),
	      progress_("CorrectionVectorSkeleton")
	{}

	void calcDynVectors(const VectorWithOffsetType& tv0,
	                    VectorWithOffsetType& tv1,
	                    VectorWithOffsetType& tv2)
	{
		const VectorWithOffsetType& phi = tv0;
		tv1 = tv2 = phi;

		VectorMatrixFieldType V(phi.sectors());
		VectorMatrixFieldType T(phi.sectors());

		VectorSizeType steps(phi.sectors());

		triDiag(phi,T,V,steps);

		VectorVectorRealType eigs(phi.sectors());

		for (SizeType ii = 0;ii < phi.sectors(); ++ii)
			PsimagLite::diag(T[ii],eigs[ii],'V');

		for (SizeType i=0;i<phi.sectors();i++) {
			VectorType sv;
			SizeType i0 = phi.sector(i);
			tv0.extract(sv,i0);
			// g.s. is included separately
			// set Aq
			// tv0.setDataInSector(sv,i0);
			// set xi
			SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
			VectorType xi(sv.size(),0),xr(sv.size(),0);

			if (tstStruct_.algorithm() == TargetParamsType::KRYLOV) {
				computeXiAndXrKrylov(xi,xr,phi,i0,V[i],T[i],eigs[i],steps[i]);
			} else {
				computeXiAndXrIndirect(xi,xr,sv,p);
			}

			tv1.setDataInSector(xi,i0);
			//set xr
			tv2.setDataInSector(xr,i0);
			DenseMatrixType V;
			getLanczosVectors(V,sv,p);
		}

		weightForContinuedFraction_ = PsimagLite::real(phi*phi);
	}

	void calcDynVectors(const VectorWithOffsetType& tv0,
	                    const VectorWithOffsetType& tv1,
	                    VectorWithOffsetType& tv2,
	                    VectorWithOffsetType& tv3)
	{
		throw PsimagLite::RuntimeError("calcDynVectors unimplemented\n");
	}

	template<typename SomeTargetingCommonType>
	void save(const SomeTargetingCommonType& targetingCommon,
	          const VectorSizeType& block,
	          PsimagLite::IoSimple::Out& io) const
	{
		if (block.size()!=1) {
			PsimagLite::String str("TargetingCorrectionVector ");
			str += "only supports blocks of size 1\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType type = tstStruct_.type();
		int fermionSign = targetingCommon.findFermionSignOfTheOperators();
		int s = (type&1) ? -1 : 1;
		int s2 = (type>1) ? -1 : 1;
		int s3 = (type&1) ? -fermionSign : 1;

		typename PostProcType::ParametersType params(ioIn_,"DynamicDmrg");
		params.Eg = targetingCommon.energy();
		params.weight = s2*weightForContinuedFraction_*s3;
		params.isign = s;
		if (ab_.size() == 0) {
			PsimagLite::OstringStream msg;
			msg<<"WARNING:  Trying to save a tridiagonal matrix with size zero.\n";
			msg<<"\tHINT: Maybe the dyn vectors were never calculated.\n";
			msg<<"\tHINT: Maybe TSPLoops is too large";
			if (params.weight != 0)
				msg<<"\n\tExpect a throw anytime now...";
			progress_.printline(msg,std::cerr);
		}

		PostProcType cf(ab_,reortho_,params);

		targetingCommon.save(block,io,cf,targetingCommon.targetVectors());
		targetingCommon.psi().save(io,"PSI");
	}

private:

	void getLanczosVectors(DenseMatrixType& V,
	                       const VectorType& sv,
	                       SizeType p)
	{
		SizeType threadId = 0;
		RealType fakeTime = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,fakeTime,threadId);
		typedef typename LanczosSolverType::LanczosMatrixType
		        LanczosMatrixType;
		LanczosMatrixType h(&model_,&modelHelper);
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
		RealType fakeTime = 0;
		typename ModelType::ModelHelperType modelHelper(p,lrs_,fakeTime,threadId);
		LanczosMatrixType h(&model_,&modelHelper);
		RealType E0 = energy_;
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
		CalcRType what(tstStruct_,energy_,eigs);

		calcR(r,what.imag(),T,V,phi,eigs,n2,i0);

		psimag::BLAS::GEMV('N',n2,n2,zone,&(T(0,0)),n2,&(r[0]),1,zzero,&(tmp[0]),1);

		xi.resize(n);
		psimag::BLAS::GEMV('N',n,n2,zone,&(V(0,0)),n,&(tmp[0]),1,zzero,&(xi[0]),1);

		calcR(r,what.real(),T,V,phi,eigs,n2,i0);

		psimag::BLAS::GEMV('N',n2,n2,zone,&(T(0,0)),n2,&(r[0]),1,zzero,&(tmp[0]),1);

		xr.resize(n);
		psimag::BLAS::GEMV('N',n,n2,zone,&(V(0,0)),n,&(tmp[0]),1,zzero,&(xr[0]),1);
	}

	void calcR(TargetVectorType& r,
	           const typename CalcRType::ActionType& whatRorI,
	           const MatrixComplexOrRealType& T,
	           const MatrixComplexOrRealType& V,
	           const VectorWithOffsetType& phi,
	           const VectorRealType&,
	           SizeType n2,
	           SizeType i0)
	{
		for (SizeType k=0;k<n2;k++) {
			ComplexOrRealType sum = 0.0;
			for (SizeType kprime=0;kprime<n2;kprime++) {
				ComplexOrRealType tmpV = calcVTimesPhi(kprime,V,phi,i0);
				sum += PsimagLite::conj(T(kprime,k))*tmpV;
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
			ret += PsimagLite::conj(V(j,kprime))*phi.fastAccess(i0,j);
		return ret;
	}

	void triDiag(const VectorWithOffsetType& phi,
	             VectorMatrixFieldType& T,
	             VectorMatrixFieldType& V,
	             typename PsimagLite::Vector<SizeType>::Type& steps)
	{
		RealType fakeTime = 0;
		typedef PsimagLite::NoPthreads<ParallelTriDiagType> ParallelizerType;
		ParallelizerType threadedTriDiag(1,0);

		ParallelTriDiagType helperTriDiag(phi,
		                                  T,
		                                  V,
		                                  steps,
		                                  lrs_,
		                                  fakeTime,
		                                  model_,
		                                  ioIn_);

		threadedTriDiag.loopCreate(phi.sectors(),helperTriDiag);
	}

	RealType dynWeightOf(VectorType& v,const VectorType& w) const
	{
		RealType sum = 0;
		for (SizeType i=0;i<v.size();i++) {
			RealType tmp = PsimagLite::real(v[i]*w[i]);
			sum += tmp*tmp;
		}
		return sum;
	}

	InputValidatorType& ioIn_;
	const TargetParamsType& tstStruct_;
	const ModelType& model_;
	const LeftRightSuperType& lrs_;
	RealType energy_;
	PsimagLite::ProgressIndicator progress_;
	TridiagonalMatrixType ab_;
	DenseMatrixRealType reortho_;
	RealType weightForContinuedFraction_;
	typename LanczosSolverType::ParametersSolverType paramsForSolver_;
}; // class CorrectionVectorSkeleton

} // namespace
/*@}*/
#endif // CORRECTION_VECTOR_SKEL_H

