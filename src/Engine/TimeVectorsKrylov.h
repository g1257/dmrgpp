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

namespace Dmrg {

template<typename TargettingParamsType,
		 typename ModelType,
		 typename WaveFunctionTransfType,
		 typename LanczosSolverType,
		 typename VectorWithOffsetType>
class TimeVectorsKrylov : public  TimeVectorsBase<
		TargettingParamsType,
		ModelType,
		WaveFunctionTransfType,
		LanczosSolverType,
		VectorWithOffsetType> {
	typedef TimeVectorsBase<TargettingParamsType,ModelType,WaveFunctionTransfType,LanczosSolverType,VectorWithOffsetType> BaseType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargettingParamsType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type TargetVectorType;
	typedef typename LanczosSolverType::TridiagonalMatrixType TridiagonalMatrixType;

public:

	TimeVectorsKrylov(RealType& currentTime,
					  const TargettingParamsType& tstStruct,
					  const VectorRealType& times,
					  typename PsimagLite::Vector<VectorWithOffsetType>::Type& targetVectors,
					  const ModelType& model,
					  const WaveFunctionTransfType& wft,
					  const LeftRightSuperType& lrs,
					  const RealType& E0)
		: currentTime_(currentTime),
		  tstStruct_(tstStruct),
		  times_(times),
		  targetVectors_(targetVectors),
		  model_(model),
		  wft_(wft),
		  lrs_(lrs),
		  E0_(E0)
	{}

	virtual void calcTimeVectors(const PairType& startEnd,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             size_t systemOrEnviron,
	                             bool allOperatorsApplied)
	{
		if (currentTime_==0 && tstStruct_.noOperator) {
			for (size_t i=0;i<times_.size();i++)
				targetVectors_[i]=phi;
			return;
		}

		calcTimeVectorsKrylov1(startEnd,Eg,phi,systemOrEnviron);
	}

private:

	void calcTimeVectorsKrylov1(const PairType& startEnd,
	                            RealType Eg,
								const VectorWithOffsetType& phi,
								size_t systemOrEnviron)
	{
		typename PsimagLite::Vector<MatrixComplexOrRealType>::Type V(phi.sectors());
		typename PsimagLite::Vector<MatrixComplexOrRealType>::Type T(phi.sectors());

		typename PsimagLite::Vector<size_t>::Type steps(phi.sectors());

		triDiag(phi,T,V,steps);

		typename PsimagLite::Vector<typename PsimagLite::Vector<RealType>::Type>::Type eigs(phi.sectors());

		for (size_t ii=0;ii<phi.sectors();ii++)
			PsimagLite::diag(T[ii],eigs[ii],'V');

		calcTargetVectors(startEnd,phi,T,V,Eg,eigs,steps,systemOrEnviron);

		//checkNorms();
	}

	//! Do not normalize states here, it leads to wrong results (!)
	void calcTargetVectors(const PairType& startEnd,
	                       const VectorWithOffsetType& phi,
						   const typename PsimagLite::Vector<MatrixComplexOrRealType>::Type& T,
						   const typename PsimagLite::Vector<MatrixComplexOrRealType>::Type& V,
						   RealType Eg,
						   const typename PsimagLite::Vector<VectorRealType>::Type& eigs,
						   typename PsimagLite::Vector<size_t>::Type steps,
						   size_t systemOrEnviron)
	{
		targetVectors_[0] = phi;
		for (size_t i=startEnd.first+1;i<startEnd.second;i++) {
			assert(i<targetVectors_.size());
			targetVectors_[i] = phi;
			// Only time differences here (i.e. times_[i] not times_[i]+currentTime_)
			calcTargetVector(targetVectors_[i],phi,T,V,Eg,eigs,i,steps);
		}
	}

	void calcTargetVector(VectorWithOffsetType& v,
	                      const VectorWithOffsetType& phi,
	                      const typename PsimagLite::Vector<MatrixComplexOrRealType>::Type& T,
	                      const typename PsimagLite::Vector<MatrixComplexOrRealType>::Type& V,
	                      RealType Eg,
	                      const typename PsimagLite::Vector<VectorRealType>::Type& eigs,
	                      size_t timeIndex,
	                      typename PsimagLite::Vector<size_t>::Type steps)
	{
		v = phi;
		for (size_t ii=0;ii<phi.sectors();ii++) {
			size_t i0 = phi.sector(ii);
			TargetVectorType r;
			calcTargetVector(r,phi,T[ii],V[ii],Eg,eigs[ii],timeIndex,steps[ii],i0);
			//std::cerr<<"TARGET FOR t="<<t<<" "<<PsimagLite::norm(r)<<" "<<norm(phi)<<"\n";
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
			size_t timeIndex,
			size_t steps,
			size_t i0)
	{
		size_t n2 = steps;
		size_t n = V.n_row();
		if (T.n_col()!=T.n_row()) throw std::runtime_error("T is not square\n");
		if (V.n_col()!=T.n_col()) throw std::runtime_error("V is not nxn2\n");
		// for (size_t j=0;j<v.size();j++) v[j] = 0; <-- harmful if v is sparse
		ComplexOrRealType zone = 1.0;
		ComplexOrRealType zzero = 0.0;

		//check1(phi,i0);
		//check2(T,V,phi,n2,i0);
		TargetVectorType tmp(n2);
		r.resize(n2);
		calcR(r,T,V,phi,Eg,eigs,timeIndex,steps,i0);
		//				std::cerr<<"TARGET FOR t="<<t<<" after calcR norm="<<PsimagLite::norm(r)<<"\n";
		psimag::BLAS::GEMV('N', n2, n2, zone, &(T(0,0)), n2, &(r[0]), 1, zzero, &(tmp[0]), 1 );
		//				std::cerr<<"TARGET FOR t="<<t<<" after S^\\dagger norm="<<PsimagLite::norm(tmp)<<"\n";
		r.resize(n);
		psimag::BLAS::GEMV('N', n,  n2, zone, &(V(0,0)), n, &(tmp[0]),1, zzero, &(r[0]),   1 );
	}

	void calcR(TargetVectorType& r,
			   const MatrixComplexOrRealType& T,
			   const MatrixComplexOrRealType& V,
			   const VectorWithOffsetType& phi,
			   RealType Eg,
			   const VectorRealType& eigs,
			   size_t timeIndex,
			   size_t n2,
			   size_t i0)
	{
		for (size_t k=0;k<n2;k++) {
			ComplexOrRealType sum = 0.0;
			for (size_t kprime=0;kprime<n2;kprime++) {
				ComplexOrRealType tmpV = calcVTimesPhi(kprime,V,phi,i0);
				sum += std::conj(T(kprime,k))*tmpV;
			}
			RealType tmp = (eigs[k]-E0_)*times_[timeIndex];
			ComplexOrRealType c = 0.0;
			PsimagLite::expComplexOrReal(c,-tmp);
			r[k] = sum * c;
		}
	}

	ComplexOrRealType calcVTimesPhi(size_t kprime,const MatrixComplexOrRealType& V,const VectorWithOffsetType& phi,
							  size_t i0) const
	{
		ComplexOrRealType ret = 0;
		size_t total = phi.effectiveSize(i0);

		for (size_t j=0;j<total;j++)
			ret += std::conj(V(j,kprime))*phi.fastAccess(i0,j);
		return ret;
	}

	void triDiag(
			const VectorWithOffsetType& phi,
			typename PsimagLite::Vector<MatrixComplexOrRealType>::Type& T,
			typename PsimagLite::Vector<MatrixComplexOrRealType>::Type& V,
			typename PsimagLite::Vector<size_t>::Type& steps)
	{
		for (size_t ii=0;ii<phi.sectors();ii++) {
			size_t i = phi.sector(ii);
			steps[ii] = triDiag(phi,T[ii],V[ii],i);
		}
	}

	size_t triDiag(const VectorWithOffsetType& phi,MatrixComplexOrRealType& T,MatrixComplexOrRealType& V,size_t i0)
	{
		size_t p = lrs_.super().findPartitionNumber(phi.offset(i0));
		typename ModelType::ModelHelperType modelHelper(p,lrs_);
				//,useReflection_);
		typename LanczosSolverType::LanczosMatrixType lanczosHelper(&model_,&modelHelper);

		typename LanczosSolverType::ParametersSolverType params;
		params.steps = model_.params().lanczosSteps;
		params.tolerance = model_.params().lanczosEps;
		params.stepsForEnergyConvergence =ProgramGlobals::MaxLanczosSteps;

		LanczosSolverType lanczosSolver(lanczosHelper,params,&V);

		TridiagonalMatrixType ab;
		size_t total = phi.effectiveSize(i0);
		TargetVectorType phi2(total);
		phi.extract(phi2,i0);
		/* std::ostringstream msg;
		msg<<"Calling tridiagonalDecomposition...\n";
		progress_.printline(msg,std::cerr);*/
		lanczosSolver.decomposition(phi2,ab);
		lanczosSolver.buildDenseMatrix(T,ab);
		//check1(V,phi2);
		return lanczosSolver.steps();
	}

	RealType& currentTime_;
	const TargettingParamsType& tstStruct_;
	const VectorRealType& times_;
	typename PsimagLite::Vector<VectorWithOffsetType>::Type& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
	const RealType& E0_;
}; //class TimeVectorsKrylov
} // namespace Dmrg
/*@}*/
#endif
