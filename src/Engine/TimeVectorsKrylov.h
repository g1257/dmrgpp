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
#include "NoPthreads.h"
#include "Parallelizer.h"

namespace Dmrg {

template<typename TargetParamsType,
         typename ModelType,
         typename WaveFunctionTransfType,
         typename LanczosSolverType,
         typename VectorWithOffsetType>
class TimeVectorsKrylov : public  TimeVectorsBase<
        TargetParamsType,
        ModelType,
        WaveFunctionTransfType,
        LanczosSolverType,
        VectorWithOffsetType> {

	typedef TimeVectorsBase<TargetParamsType,
	ModelType,
	WaveFunctionTransfType,
	LanczosSolverType,
	VectorWithOffsetType> BaseType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType
	BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef ParallelTriDiag<ModelType,LanczosSolverType,VectorWithOffsetType>
	ParallelTriDiagType;
	typedef typename ParallelTriDiagType::MatrixComplexOrRealType
	MatrixComplexOrRealType;
	typedef typename ParallelTriDiagType::TargetVectorType TargetVectorType;
	typedef typename ParallelTriDiagType::VectorMatrixFieldType
	VectorMatrixFieldType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;
	typedef typename ModelType::InputValidatorType InputValidatorType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type
	VectorVectorWithOffsetType;
	typedef typename PsimagLite::Vector<VectorRealType>::Type VectorVectorRealType;

public:

	TimeVectorsKrylov(const RealType& currentTime,
	                  const TargetParamsType& tstStruct,
	                  const VectorRealType& times,
	                  VectorVectorWithOffsetType& targetVectors,
	                  const ModelType& model,
	                  const WaveFunctionTransfType& wft,
	                  const LeftRightSuperType& lrs,
	                  const RealType& E0,
	                  InputValidatorType& ioIn)
	    : currentTime_(currentTime),
	      tstStruct_(tstStruct),
	      times_(times),
	      targetVectors_(targetVectors),
	      model_(model),
	      wft_(wft),
	      lrs_(lrs),
	      E0_(E0),
	      ioIn_(ioIn),
	      timeHasAdvanced_(true)
	{}

	virtual void calcTimeVectors(const PairType& startEnd,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             SizeType systemOrEnviron,
	                             bool,
	                             const PsimagLite::Vector<SizeType>::Type&)
	{
		if (currentTime_==0 && tstStruct_.noOperator() && tstStruct_.skipTimeZero()) {
			for (SizeType i=0;i<times_.size();i++)
				targetVectors_[i]=phi;
			return;
		}

		bool returnEarly = true;

		for (SizeType i = 0; i < targetVectors_.size(); ++i) {
			if (targetVectors_[i].size() == 0) {
				returnEarly = false;
				break;
			}
		}

		assert(0 < targetVectors_.size());
		targetVectors_[0] = phi;
		if (timeHasAdvanced_) returnEarly = false;
		if (returnEarly && ProgramGlobals::TST_FAST) return;

		VectorMatrixFieldType V(phi.sectors());
		VectorMatrixFieldType T(phi.sectors());

		typename PsimagLite::Vector<SizeType>::Type steps(phi.sectors());

		triDiag(phi,T,V,steps);

		VectorVectorRealType eigs(phi.sectors());

		for (SizeType ii=0;ii<phi.sectors();ii++)
			PsimagLite::diag(T[ii],eigs[ii],'V');

		calcTargetVectors(startEnd,phi,T,V,Eg,eigs,steps,systemOrEnviron);

		//checkNorms();
		timeHasAdvanced_ = false;
	}

	void timeHasAdvanced()
	{
		timeHasAdvanced_ = true;
	}

private:

	//! Do not normalize states here, it leads to wrong results (!)
	void calcTargetVectors(const PairType& startEnd,
	                       const VectorWithOffsetType& phi,
	                       const VectorMatrixFieldType& T,
	                       const VectorMatrixFieldType& V,
	                       RealType Eg,
	                       const VectorVectorRealType& eigs,
	                       typename PsimagLite::Vector<SizeType>::Type steps,
	                       SizeType)
	{
		for (SizeType i=startEnd.first+1;i<startEnd.second;i++) {
			assert(i<targetVectors_.size());
			targetVectors_[i] = phi;
			// Only time differences here (i.e. times_[i] not times_[i]+currentTime_)
			calcTargetVector(targetVectors_[i],phi,T,V,Eg,eigs,i,steps);
		}
	}

	void calcTargetVector(VectorWithOffsetType& v,
	                      const VectorWithOffsetType& phi,
	                      const VectorMatrixFieldType& T,
	                      const VectorMatrixFieldType& V,
	                      RealType Eg,
	                      const VectorVectorRealType& eigs,
	                      SizeType timeIndex,
	                      typename PsimagLite::Vector<SizeType>::Type steps)
	{
		v = phi;
		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i0 = phi.sector(ii);
			TargetVectorType r;
			calcTargetVector(r,phi,T[ii],V[ii],Eg,eigs[ii],timeIndex,steps[ii],i0);
			v.setDataInSector(r,i0);
		}
	}

	void calcTargetVector(
	        TargetVectorType& r,
	        const VectorWithOffsetType& phi,
	        const MatrixComplexOrRealType& T,
	        const MatrixComplexOrRealType& V,
	        RealType Eg,
	        const VectorRealType& eigs,
	        SizeType timeIndex,
	        SizeType steps,
	        SizeType i0)
	{
		SizeType n2 = steps;
		SizeType n = V.n_row();
		if (T.n_col()!=T.n_row()) throw PsimagLite::RuntimeError("T is not square\n");
		if (V.n_col()!=T.n_col()) throw PsimagLite::RuntimeError("V is not nxn2\n");
		// for (SizeType j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
		ComplexOrRealType zone = 1.0;
		ComplexOrRealType zzero = 0.0;

		//check1(phi,i0);
		//check2(T,V,phi,n2,i0);
		TargetVectorType tmp(n2);
		r.resize(n2);
		calcR(r,T,V,phi,Eg,eigs,timeIndex,steps,i0);
		psimag::BLAS::GEMV('N',n2,n2,zone,&(T(0,0)),n2,&(r[0]),1,zzero,&(tmp[0]),1);
		r.resize(n);
		psimag::BLAS::GEMV('N',n,n2,zone,&(V(0,0)),n,&(tmp[0]),1,zzero,&(r[0]),1);
	}

	void calcR(TargetVectorType& r,
	           const MatrixComplexOrRealType& T,
	           const MatrixComplexOrRealType& V,
	           const VectorWithOffsetType& phi,
	           RealType,
	           const VectorRealType& eigs,
	           SizeType timeIndex,
	           SizeType n2,
	           SizeType i0)
	{
		RealType timeDirection = tstStruct_.timeDirection();

		for (SizeType k=0;k<n2;k++) {
			ComplexOrRealType sum = 0.0;
			for (SizeType kprime=0;kprime<n2;kprime++) {
				ComplexOrRealType tmpV = calcVTimesPhi(kprime,V,phi,i0);
				sum += std::conj(T(kprime,k))*tmpV;
			}

			RealType tmp = (eigs[k]-E0_)*times_[timeIndex]*timeDirection;
			ComplexOrRealType c = 0.0;
			PsimagLite::expComplexOrReal(c,-tmp);
			r[k] = sum * c;
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
		PsimagLite::String options = model_.params().options;
		bool cTridiag = (options.find("concurrenttridiag") != PsimagLite::String::npos);

		if (cTridiag) {
			typedef PsimagLite::Parallelizer<ParallelTriDiagType> ParallelizerType;
			ParallelizerType threadedTriDiag(PsimagLite::Concurrency::npthreads,
			                                 PsimagLite::MPI::COMM_WORLD);

			ParallelTriDiagType helperTriDiag(phi,T,V,steps,lrs_,currentTime_,model_,ioIn_);

			threadedTriDiag.loopCreate(phi.sectors(),helperTriDiag);
		} else {
			typedef PsimagLite::NoPthreads<ParallelTriDiagType> ParallelizerType;
			ParallelizerType threadedTriDiag(1,0);

			ParallelTriDiagType helperTriDiag(phi,T,V,steps,lrs_,currentTime_,model_,ioIn_);

			threadedTriDiag.loopCreate(phi.sectors(),helperTriDiag);
		}
	}

	const RealType& currentTime_;
	const TargetParamsType& tstStruct_;
	const VectorRealType& times_;
	VectorVectorWithOffsetType& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
	const RealType& E0_;
	InputValidatorType& ioIn_;
	bool timeHasAdvanced_;
}; //class TimeVectorsKrylov
} // namespace Dmrg
/*@}*/
#endif

