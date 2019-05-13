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

/*! \file TimeVectorsKrylov.h
 *
 *
 */

#ifndef TIME_VECTORS_KRYLOV
#define TIME_VECTORS_KRYLOV
#include <iostream>
#include <vector>
#include "TimeVectorsBase.h"
#include "ParallelTriDiag.h"
#include "NoPthreadsNg.h"
#include "Parallelizer.h"
#include "KrylovHelper.h"

namespace Dmrg {

template<typename TargetParamsType,
         typename ModelType,
         typename WaveFunctionTransfType,
         typename LanczosSolverType,
         typename VectorWithOffsetType_>
class TimeVectorsKrylov : public  TimeVectorsBase<
        TargetParamsType,
        ModelType,
        WaveFunctionTransfType,
        LanczosSolverType,
        VectorWithOffsetType_> {

public:

	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef TimeVectorsBase<TargetParamsType,
	ModelType,
	WaveFunctionTransfType,
	LanczosSolverType,
	VectorWithOffsetType> BaseType;
	typedef TimeVectorsKrylov<TargetParamsType,
	ModelType,
	WaveFunctionTransfType,
	LanczosSolverType,
	VectorWithOffsetType> ThisType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef ParallelTriDiag<ModelType,LanczosSolverType,VectorWithOffsetType> ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::TargetVectorType VectorType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType VectorMatrixFieldType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;

	struct Action {

		typedef typename ThisType::MatrixComplexOrRealType MatrixComplexOrRealType;
		typedef typename ThisType::VectorWithOffsetType VectorWithOffsetType;
		typedef typename ThisType::VectorRealType VectorRealType;
		typedef typename ModelType::SolverParamsType SolverParamsType;
	};

	typedef KrylovHelper<Action> KrylovHelperType;

	TimeVectorsKrylov(const SizeType& currentTimeStep,
	                  const TargetParamsType& tstStruct,
	                  const VectorRealType& times,
	                  VectorVectorWithOffsetType& targetVectors,
	                  const ModelType& model,
	                  const WaveFunctionTransfType& wft,
	                  const LeftRightSuperType& lrs,
	                  InputValidatorType& ioIn)
	    : currentTimeStep_(currentTimeStep),
	      tstStruct_(tstStruct),
	      times_(times),
	      targetVectors_(targetVectors),
	      model_(model),
	      wft_(wft),
	      lrs_(lrs),
	      ioIn_(ioIn),
	      timeHasAdvanced_(false),
	      krylovHelper_(model.params())
	{}

	virtual void calcTimeVectors(const PsimagLite::Vector<SizeType>::Type& indices,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             typename BaseType::ExtraData* = 0)
	{
		const SizeType n = indices.size();
		if (currentTimeStep_ == 0 && tstStruct_.noOperator() && tstStruct_.skipTimeZero()) {
			for (SizeType i = 0; i < n; ++i) {
				const SizeType ii = indices[i];
				targetVectors_[ii] = phi;
			}
		}

		const VectorWithOffsetType* ptr0 = &(targetVectors_[indices[0]]);
		const VectorWithOffsetType* ptr1 = &phi;
		if (ptr0 != ptr1)
			targetVectors_[indices[0]] = phi;

		if (times_.size() == 1 && fabs(times_[0])<1e-10) return;

		VectorMatrixFieldType V(phi.sectors());
		VectorMatrixFieldType T(phi.sectors());

		typename PsimagLite::Vector<SizeType>::Type steps(phi.sectors());

		triDiag(phi,T,V,steps);

		VectorVectorRealType eigs(phi.sectors());

		for (SizeType ii=0;ii<phi.sectors();ii++)
			PsimagLite::diag(T[ii],eigs[ii],'V');

		calcTargetVectors(indices, phi, T, V, Eg, eigs, steps);

		//checkNorms();
		timeHasAdvanced_ = false;
	}

	void timeHasAdvanced()
	{
		timeHasAdvanced_ = true;
	}

	RealType time() const
	{
		return currentTimeStep_*tstStruct_.tau();
	}

private:

	//! Do not normalize states here, it leads to wrong results (!)
	void calcTargetVectors(typename PsimagLite::Vector<SizeType>::Type indices,
	                       const VectorWithOffsetType& phi,
	                       const VectorMatrixFieldType& T,
	                       const VectorMatrixFieldType& V,
	                       RealType Eg,
	                       const VectorVectorRealType& eigs,
	                       typename PsimagLite::Vector<SizeType>::Type steps)
	{
		for (SizeType i = 1; i < indices.size(); ++i) {
			const SizeType ii = indices[i];
			assert(ii < targetVectors_.size());
			targetVectors_[ii] = phi;
			// Only time differences here (i.e. times_[i] not times_[i]+currentTime_)
			calcTargetVector(targetVectors_[ii], phi, T, V, Eg, eigs, steps, i);
		}
	}

	void calcTargetVector(VectorWithOffsetType& v,
	                      const VectorWithOffsetType& phi,
	                      const VectorMatrixFieldType& T,
	                      const VectorMatrixFieldType& V,
	                      RealType Eg,
	                      const VectorVectorRealType& eigs,
	                      typename PsimagLite::Vector<SizeType>::Type steps,
	                      SizeType timeIndex)
	{
		v = phi;
		for (SizeType ii = 0;ii < phi.sectors(); ++ii) {
			const RealType time = times_[timeIndex];
			const RealType timeDirection = tstStruct_.timeDirection();
			const VectorRealType& eigsii = eigs[ii];
			auto action = [eigsii, Eg, time, timeDirection](SizeType k)
			{
				RealType tmp = (eigsii[k]-Eg)*time*timeDirection;
				ComplexOrRealType c = 0.0;
				PsimagLite::expComplexOrReal(c, -tmp);
				return c;
			};

			SizeType i0 = phi.sector(ii);
			VectorType r;
			calcTargetVector(r, phi, T[ii], V[ii], action, steps[ii], i0);
			v.setDataInSector(r, i0);
		}
	}

	template<typename SomeLambdaType>
	void calcTargetVector(VectorType& r,
	                      const VectorWithOffsetType& phi,
	                      const MatrixComplexOrRealType& T,
	                      const MatrixComplexOrRealType& V,
	                      const SomeLambdaType& action,
	                      SizeType steps,
	                      SizeType i0)
	{
		SizeType n2 = steps;
		SizeType n = V.rows();
		if (T.cols()!=T.rows()) throw PsimagLite::RuntimeError("T is not square\n");
		if (V.cols()!=T.cols()) throw PsimagLite::RuntimeError("V is not nxn2\n");
		// for (SizeType j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
		ComplexOrRealType zone = 1.0;
		ComplexOrRealType zzero = 0.0;

		//check1(phi,i0);
		//check2(T,V,phi,n2,i0);
		VectorType tmp(n2);
		r.resize(n2);
		krylovHelper_.calcR(r, action, T, V, phi, steps, i0);
		psimag::BLAS::GEMV('N',n2,n2,zone,&(T(0,0)),n2,&(r[0]),1,zzero,&(tmp[0]),1);
		r.resize(n);
		psimag::BLAS::GEMV('N',n,n2,zone,&(V(0,0)),n,&(tmp[0]),1,zzero,&(r[0]),1);
	}

	void triDiag(const VectorWithOffsetType& phi,
	             VectorMatrixFieldType& T,
	             VectorMatrixFieldType& V,
	             typename PsimagLite::Vector<SizeType>::Type& steps)
	{
		typedef PsimagLite::NoPthreadsNg<ParallelTriDiagType> ParallelizerType;
		ParallelizerType threadedTriDiag(PsimagLite::CodeSectionParams(1));

		ParallelTriDiagType helperTriDiag(phi,T,V,steps,lrs_,time(),model_,ioIn_);

		threadedTriDiag.loopCreate(helperTriDiag);
	}

	const SizeType& currentTimeStep_;
	const TargetParamsType& tstStruct_;
	const VectorRealType& times_;
	VectorVectorWithOffsetType& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
	InputValidatorType& ioIn_;
	bool timeHasAdvanced_;
	KrylovHelperType krylovHelper_;
}; //class TimeVectorsKrylov
} // namespace Dmrg
/*@}*/
#endif

