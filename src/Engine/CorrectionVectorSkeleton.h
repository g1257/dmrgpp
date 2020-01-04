/*
Copyright (c) 2009-2016-2018-2019, UT-Battelle, LLC
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
#include "FreqEnum.h"
#include "NoPthreadsNg.h"
#include "TridiagRixsStatic.h"
#include "KrylovHelper.h"

namespace Dmrg {

template<typename LanczosSolverType_,
         typename VectorWithOffsetType_,
         typename TargetingBaseType,
         typename TargetParamsType>
class CorrectionVectorSkeleton {

public:

	typedef CorrectionVectorSkeleton<LanczosSolverType_,
	VectorWithOffsetType_,
	TargetingBaseType,
	TargetParamsType> ThisType;
	typedef LanczosSolverType_ LanczosSolverType;
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
	typedef typename TargetingBaseType::TargetingCommonType TargetingCommonType;
	typedef typename TargetingCommonType::TimeSerializerType TimeSerializerType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef PsimagLite::Matrix<typename VectorType::value_type> DenseMatrixType;
	typedef PsimagLite::Matrix<RealType> DenseMatrixRealType;
	typedef typename LanczosSolverType::PostProcType PostProcType;
	typedef typename LanczosSolverType::MatrixType LanczosMatrixType;
	typedef CorrectionVectorFunction<LanczosMatrixType,
	TargetParamsType> CorrectionVectorFunctionType;
	typedef ParallelTriDiag<ModelType, LanczosSolverType, VectorWithOffsetType>
	ParallelTriDiagType;
	typedef TridiagRixsStatic<ModelType, LanczosSolverType, VectorWithOffsetType>
	TridiagRixsStaticType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;
	typedef typename ModelType::InputValidatorType InputValidatorType;

	class CalcR {

	public:

		class Action {

		public:

			typedef typename ThisType::ModelType::SolverParamsType SolverParamsType;
			typedef typename ThisType::MatrixComplexOrRealType MatrixComplexOrRealType;
			typedef typename ThisType::VectorWithOffsetType VectorWithOffsetType;
			typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

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

	typedef KrylovHelper<typename CalcR::ActionType> KrylovHelperType;

	CorrectionVectorSkeleton(InputValidatorType& ioIn,
	                         const TargetParamsType& tstStruct,
	                         const ModelType& model,
	                         const LeftRightSuperType& lrs,
	                         const RealType& energy)
	    : ioIn_(ioIn),
	      tstStruct_(tstStruct),
	      model_(model),
	      lrs_(lrs),
	      energy_(energy),
	      progress_("CorrectionVectorSkeleton"),
	      krylovHelper_(model.params())
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
			//tv0.setDataInSector(sv,i0);
			// set xi
			SizeType p = lrs_.super().findPartitionNumber(phi.offset(i0));
			VectorType xi(sv.size(),0),xr(sv.size(),0);

			if (tstStruct_.algorithm() == TargetParamsType::BaseType::AlgorithmEnum::KRYLOV) {
				computeXiAndXrKrylov(xi,xr,phi,i0,V[i],T[i],eigs[i],steps[i]);
			} else {
				computeXiAndXrIndirect(xi,xr,sv,p);
			}

			tv1.setDataInSector(xi,i0);
			//set xr
			tv2.setDataInSector(xr,i0);
		}

		weightForContinuedFraction_ = PsimagLite::real(phi*phi);
	}

	void calcDynVectors(const VectorWithOffsetType& tv0,
	                    const VectorWithOffsetType& tv1,
	                    VectorWithOffsetType& tv2,
	                    VectorWithOffsetType& tv3)
	{
		VectorWithOffsetType tv4;
		calcDynVectors(tv0,tv4,tv2);
		VectorWithOffsetType tv5;
		calcDynVectors(tv1,tv5,tv3);
		tv2 += tv5;
		tv3 += (-1.0)*tv4;
	}

private:

	void computeXiAndXrIndirect(VectorType& xi,
	                            VectorType& xr,
	                            const VectorType& sv,
	                            SizeType p)
	{
		if (tstStruct_.omega().first != PsimagLite::FREQ_REAL)
			throw PsimagLite::RuntimeError("Matsubara only with KRYLOV\n");

		const RealType fakeTime = 0;
		typename ModelHelperType::Aux aux(p, lrs_);
		typename ModelType::HamiltonianConnectionType hc(lrs_,
		                                                 ModelType::modelLinks(),
		                                                 fakeTime,
		                                                 model_.superOpHelper());
		LanczosMatrixType h(model_, hc, aux);
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

		ComplexOrRealType zone = 1.0;
		ComplexOrRealType zzero = 0.0;

		TargetVectorType tmp(n2);
		VectorType r(n2);
		CalcR what(tstStruct_, energy_, eigs);

		krylovHelper_.calcR(r, what.imag(), T, V, phi, n2, i0);

		psimag::BLAS::GEMV('N',n2,n2,zone,&(T(0,0)),n2,&(r[0]),1,zzero,&(tmp[0]),1);

		xi.resize(n);
		psimag::BLAS::GEMV('N',n,n2,zone,&(V(0,0)),n,&(tmp[0]),1,zzero,&(xi[0]),1);

		krylovHelper_.calcR(r, what.real(), T, V, phi, n2, i0);

		psimag::BLAS::GEMV('N',n2,n2,zone,&(T(0,0)),n2,&(r[0]),1,zzero,&(tmp[0]),1);

		xr.resize(n);
		psimag::BLAS::GEMV('N',n,n2,zone,&(V(0,0)),n,&(tmp[0]),1,zzero,&(xr[0]),1);
	}

	void triDiag(const VectorWithOffsetType& phi,
	             VectorMatrixFieldType& T,
	             VectorMatrixFieldType& V,
	             VectorSizeType& steps)
	{
		RealType fakeTime = 0;
		typedef PsimagLite::NoPthreadsNg<ParallelTriDiagType> ParallelizerType;
		ParallelizerType threadedTriDiag(PsimagLite::CodeSectionParams(1));

		ParallelTriDiagType helperTriDiag(phi,
		                                  T,
		                                  V,
		                                  steps,
		                                  lrs_,
		                                  fakeTime,
		                                  model_,
		                                  ioIn_);

		threadedTriDiag.loopCreate(helperTriDiag);
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
	const RealType& energy_;
	PsimagLite::ProgressIndicator progress_;
	RealType weightForContinuedFraction_;
	KrylovHelperType krylovHelper_;
}; // class CorrectionVectorSkeleton

} // namespace
/*@}*/
#endif // CORRECTION_VECTOR_SKEL_H

