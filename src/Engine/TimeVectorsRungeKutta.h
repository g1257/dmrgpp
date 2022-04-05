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

/*! \file TimeVectorsRungeKutta.h
 *
 *
 */

#ifndef TIME_VECTORS_RUNGE_KUTTA
#define TIME_VECTORS_RUNGE_KUTTA
#include <iostream>
#include "RungeKutta.h"
#include "TimeVectorsBase.h"

namespace Dmrg {

template<typename RealType>
RealType minusOneOrMinusI(RealType)
{
	return -1;
}

template<typename RealType>
std::complex<RealType> minusOneOrMinusI(std::complex<RealType>)
{
	return std::complex<RealType>(0.0,-1.0);
}

template<typename TargetParamsType,
		 typename ModelType,
		 typename WaveFunctionTransfType,
		 typename LanczosSolverType,
		 typename VectorWithOffsetType>
class TimeVectorsRungeKutta : public  TimeVectorsBase<
		TargetParamsType,
		ModelType,
		WaveFunctionTransfType,
		LanczosSolverType,
		VectorWithOffsetType> {

	typedef TimeVectorsBase<TargetParamsType,
	ModelType,WaveFunctionTransfType,
	LanczosSolverType,
	VectorWithOffsetType> BaseType;
	typedef typename BaseType::PairType PairType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexOrRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef VectorComplexOrRealType TargetVectorType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType*>::Type VectorVectorWithOffsetType;

public:

	TimeVectorsRungeKutta(const TargetParamsType& tstStruct,
						  VectorVectorWithOffsetType& targetVectors,
						  const ModelType& model,
						  const WaveFunctionTransfType& wft,
						  const LeftRightSuperType& lrs)
		: BaseType(model, lrs, wft, "rungekutta"),
	      progress_("TimeVectorsRungeKutta"),
	      tstStruct_(tstStruct),
		  targetVectors_(targetVectors),
		  model_(model),
		  wft_(wft),
		  lrs_(lrs)
	{}

	virtual void calcTimeVectors(const VectorSizeType& indices,
	                             RealType Eg,
	                             const VectorWithOffsetType& phi,
	                             const typename BaseType::ExtraData&)
	{
		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"using RungeKutta";

		RealType norma = norm(phi);
		if (norma<1e-10) return;
		msg<<" Norm of phi= "<<norma;
		progress_.printline(msgg, std::cout);

		// set non-zero sectors
		for (SizeType i=0;i<tstStruct_.times().size();i++) *targetVectors_[i] = phi;

		for (SizeType ii=0;ii<phi.sectors();ii++) {
			SizeType i = phi.sector(ii);
			calcTimeVectors(indices, Eg, phi, i);
		}
	}

private:

	class FunctionForRungeKutta {

	public:

		FunctionForRungeKutta(const RealType& E0,
		                      const RealType &timeDirection,
		                      const LeftRightSuperType& lrs,
		                      RealType currentTime,
		                      const ModelType& model,
		                      const VectorWithOffsetType& phi,
		                      SizeType i0)
			: E0_(E0),
		      timeDirection_(timeDirection),
			  p_(lrs.super().findPartitionNumber(phi.offset(i0))),
		      aux_(p_, lrs),
		      hc_(lrs, ModelType::modelLinks(), currentTime, model.superOpHelper()),
			  lanczosHelper_(model, hc_, aux_)
		{}

		TargetVectorType operator()(const RealType&,const TargetVectorType& y) const
		{
			TargetVectorType x(y.size());
			lanczosHelper_.matrixVectorProduct(x,y);
			for (SizeType i=0;i<x.size();i++) x[i] -= E0_*y[i];
			ComplexOrRealType icomplex = minusOneOrMinusI(static_cast<ComplexOrRealType>(0));
			ComplexOrRealType tmp2 = timeDirection_*icomplex;
			TargetVectorType x2;
			x2 <= tmp2*x;
			return x2;
		}

	private:

		RealType E0_;
		RealType timeDirection_;
		SizeType p_;
		typename ModelHelperType::Aux aux_;
		typename ModelType::HamiltonianConnectionType hc_;
		typename LanczosSolverType::MatrixType lanczosHelper_;
	}; // FunctionForRungeKutta

	void calcTimeVectors(const VectorSizeType& indices,
	                     RealType Eg,
	                     const VectorWithOffsetType& phi,
	                     SizeType i0)
	{
		const VectorRealType& times = tstStruct_.times();
		SizeType total = phi.effectiveSize(i0);
		TargetVectorType phi0(total);
		phi.extract(phi0,i0);
		FunctionForRungeKutta f(Eg,
		                        tstStruct_.timeDirection(),
		                        lrs_,
		                        this->time(),
		                        model_,
		                        phi,
		                        i0);

		RealType epsForRK = tstStruct_.tau()/(times.size()-1.0);
		PsimagLite::RungeKutta<RealType,FunctionForRungeKutta,TargetVectorType>
		        rungeKutta(f,epsForRK);

		typename PsimagLite::Vector<TargetVectorType>::Type result;
		rungeKutta.solve(result,0.0,times.size(),phi0);
		assert(result.size()==times.size());

		const SizeType n = indices.size();
		for (SizeType i = 0; i < n; ++i) {
			const SizeType ii = indices[i];
			assert(ii < targetVectors_.size());
			targetVectors_[ii]->setDataInSector(result[i], i0);
		}
	}

	PsimagLite::ProgressIndicator progress_;
	const TargetParamsType& tstStruct_;
	VectorVectorWithOffsetType& targetVectors_;
	const ModelType& model_;
	const WaveFunctionTransfType& wft_;
	const LeftRightSuperType& lrs_;
}; //class TimeVectorsRungeKutta
} // namespace Dmrg
/*@}*/
#endif
