
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2011, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************

*/
/** \ingroup PsimagLite */
/*@{*/

/*! \file ChebyshevSerializer.h
 *
 * This class saves the moments of a Chebyshev expansion so that
 * they can later be processed (with KernelPolynomial, for example)
 *
 * The interface follows ContinuedFraction so that it can be
 * used interchangably
 *
 */

#ifndef CHEBYSHEV_SERIALIZER_H
#define CHEBYSHEV_SERIALIZER_H
#include <iostream>
#include "TypeToString.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include <stdexcept>
#include "ParametersForSolver.h"
#include "PlotParams.h"
#include "ChebyshevFunction.h"
#include <cassert>
#include "String.h"

namespace PsimagLite {

template<typename RealType>
struct KernelPolynomialParameters {

	enum {JACKSON,LORENTZ,DIRICHLET};

	KernelPolynomialParameters(SizeType type1,
	                           SizeType cutoff1,
	                           const RealType& lambda1)
	    : type(type1),cutoff(cutoff1),lambda(lambda1)
	{}

	SizeType type;
	SizeType cutoff;
	RealType lambda;
}; // struct KernelPolynomialParameters

template<typename VectorType_>
class ChebyshevSerializer  {

	typedef typename VectorType_::value_type VectorElementType;
	typedef typename PsimagLite::Real<VectorElementType>::Type RealType;

	static const String stringMarker_;

public:

	typedef VectorType_ VectorType;
	typedef typename VectorType::value_type FieldType;
	typedef Matrix<FieldType> MatrixType;
	typedef typename Vector<std::pair<RealType,RealType> >::Type PlotDataType;
	typedef PlotParams<RealType> PlotParamsType;
	typedef ParametersForSolver<RealType> ParametersType;
	typedef KernelPolynomialParameters<RealType> KernelParametersType;

	ChebyshevSerializer(const VectorType& ab,const ParametersType& params)
	    : progress_("ChebyshevSerializer"),
	      moments_(ab),
	      params_(params)
	{}

	template<typename IoInputType>
	ChebyshevSerializer(IoInputType& io)
	    : progress_("ChebyshevSerializer")
	{
		io.readline(params_.Eg,"#ChebyshevEnergy=");
		io.readline(params_.oneOverA,"#ChebyshevOneOverA=");
		io.readline(params_.b,"#ChebyshevB");
		io.read(moments_,"#ChebyshevMoments");
	}

	template<typename IoOutputType>
	void save(IoOutputType& io) const
	{
		String s(stringMarker_);
		io.print(s);
		io.print("#ChebyshevEnergy=",params_.Eg);
		io.print("#ChebyshevOneOverA=",params_.oneOverA);
		io.print("#ChebyshevB=",params_.b);

		io.printVector(moments_,"#ChebyshevMoments");
	}

	static const String& stringMarker() { return stringMarker_; }

	void plot(PlotDataType& result,
	          const PlotParamsType& params,
	          const KernelParametersType& kernelParams) const
	{
		SizeType cutoff = kernelParams.cutoff;
		if (cutoff==0 || moments_.size()<cutoff) cutoff = moments_.size();
		typename Vector<RealType>::Type gn(cutoff,1.0);
		initKernel(gn,kernelParams);

		typename Vector<RealType>::Type gnmun(gn.size());
		computeGnMuN(gnmun,gn);

		SizeType counter = 0;
		SizeType n = SizeType((params.omega2 - params.omega1)/params.deltaOmega);
		if (result.size()==0) result.resize(n);
		RealType offset = params_.Eg;
		std::cerr<<"gn[0]="<<gn[0]<<" gn[5]="<<gn[5]<<"\n";
		for (RealType omega=params.omega1;omega<params.omega2;omega+=params.deltaOmega) {
			RealType x = (omega+offset-params_.b)*params_.oneOverA;

			RealType den = (x>1.0 || x<-1.0) ? 0.0 : sqrt(1.0 - x*x);
			RealType tmp = (fabs(den)>1e-6) ? calcF(x,gnmun)/den : 0.0;
			std::pair<RealType,RealType> p(omega,tmp);
			result[counter++] = p;

			if (counter>=result.size()) break;
		}
	}

	//! Cases:
	//! (1) < phi0|A (z+(E0-e_k))^{-1}|A^\dagger|phi0> and
	//! (2) < phi0|A^\dagger (z-(E0-e_k))^{-1}|A|phi0>
	//! (There are actually 4 cases for the off-diagonal gf because
	//! A has two cases:
	//! (1) A = c_i + c_j and
	//! (2) A = c_i - c_j
	RealType iOfOmega(const RealType& z,RealType offset,int isign) const

	{
		throw RuntimeError("iOfOmega: unimplemented\n");
	}

private:

	RealType calcF(const RealType& x,
	               const typename Vector<RealType>::Type& gnmn) const
	{
		RealType sum = 0.5*gnmn[0];
		for (SizeType i=1;i<gnmn.size();i++) sum += gnmn[i]*chebyshev_(i,x);
		return 2.0*sum;
	}

	void computeGnMuN(typename Vector<RealType>::Type& gnmn,
	                  typename Vector<RealType>::Type& gn) const
	{
		for (SizeType i=0;i<gnmn.size();i++) gnmn[i] = moments_[i] * gn[i];
	}

	void initKernel(typename Vector<RealType>::Type& gn,
	                const KernelParametersType& kernelParams) const
	{
		switch (kernelParams.type) {
		case KernelParametersType::JACKSON:
			initKernelJackson(gn);
			break;
		case KernelParametersType::LORENTZ:
			initKernelLorentz(gn,kernelParams.lambda);
			break;
		case KernelParametersType::DIRICHLET:
			break;
		default:
			assert(false);
		}
	}

	void initKernelJackson(typename Vector<RealType>::Type& gn) const
	{
		SizeType nPlus1 = gn.size()+1;
		RealType cot1 = 1.0/tan(M_PI/nPlus1);
		for (SizeType i=0;i<gn.size();i++) {
			gn[i] = (nPlus1-i)*cos(M_PI*i/nPlus1)+sin(M_PI*i/nPlus1)*cot1;
			gn[i] /= nPlus1;
		}
	}

	void initKernelLorentz(typename Vector<RealType>::Type& gn,
	                       const RealType& lambda) const
	{
		RealType nreal = gn.size();
		RealType sinhlambda = sinh(lambda);
		for (SizeType i=0;i<gn.size();i++) {
			gn[i] = sinh(lambda*(1-i/nreal))/sinhlambda;
		}
	}

	ProgressIndicator progress_;
	typename Vector<RealType>::Type moments_;
	ParametersType params_;
	ChebyshevFunction<RealType> chebyshev_;
}; // class ChebyshevSerializer

template<typename VectorType>
const String ChebyshevSerializer<VectorType>::stringMarker_ = "#ChebyshevSerializerMarker";
} // namespace PsimagLite 
/*@}*/
#endif  //CHEBYSHEV_SERIALIZER_H
