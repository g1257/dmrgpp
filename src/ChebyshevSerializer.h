
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

namespace PsimagLite {
	
	template<typename RealType>
	struct KernelPolynomialParameters {

		enum {JACKSON,LORENTZ,DIRICHLET};

		KernelPolynomialParameters(size_t type1,size_t cutoff1,const RealType& lambda1)
		: type(type1),cutoff(cutoff1),lambda(lambda1)
		{}

		size_t type;
		size_t cutoff;
		RealType lambda;
	}; // struct KernelPolynomialParameters
	
	template<typename RealType,typename VectorType_>
	class ChebyshevSerializer  {

		static const std::string stringMarker_;

	public:

		typedef VectorType_ VectorType;
		typedef typename VectorType::value_type FieldType;
		typedef Matrix<FieldType> MatrixType;
		typedef std::vector<std::pair<RealType,RealType> > PlotDataType;
		typedef PlotParams<RealType> PlotParamsType;
		typedef ParametersForSolver<RealType> ParametersType;
		typedef KernelPolynomialParameters<RealType> KernelParametersType;

		ChebyshevSerializer(const VectorType& ab,const ParametersType& params)
		: progress_("ChebyshevSerializer",0),
		  moments_(ab),
		  params_(params)
		{}

// 		ChebyshevSerializer() 
// 		: progress_("ChebyshevSerializer",0), ab_(),Eg_(0),weight_(0),isign_(1)
// 		{ }

		template<typename IoInputType>
		ChebyshevSerializer(IoInputType& io)
		: progress_("ChebyshevSerializer",0)
		{
 			io.readline(params_.Eg,"#ChebyshevEnergy=");
			io.readline(params_.oneOverA,"#ChebyshevOneOverA=");
			io.readline(params_.b,"#ChebyshevB");
			io.read(moments_,"#ChebyshevMoments");
		}

		template<typename IoOutputType>
		void save(IoOutputType& io) const
		{
			std::string s(stringMarker_);
			io.printline(s);
			s = "#ChebyshevEnergy=";
			s += ttos(params_.Eg);
			io.printline(s);
			s = "#ChebyshevOneOverA=" + ttos(params_.oneOverA);
			io.printline(s);
			s = "#ChebyshevB=" + ttos(params_.b);
			io.printline(s);

			io.printVector(moments_,"#ChebyshevMoments");
		}
		
		static const std::string& stringMarker() { return stringMarker_; }
		

// 		void set(const VectorType& ab,
// 		         RealType Eg,
// 		         RealType weight = 0,
// 		         int isign = 0)
// 		{
// 			moments_=ab;
// 			//gatherEvenAndOdd(ab);
// 		}

		void plot(PlotDataType& result,
		          const PlotParamsType& params,
		          const KernelParametersType& kernelParams) const
		{
			size_t cutoff = kernelParams.cutoff;
			if (cutoff==0 || moments_.size()<cutoff) cutoff = moments_.size();
			std::vector<RealType> gn(cutoff,1.0);
			initKernel(gn,kernelParams);

			std::vector<RealType> gnmun(gn.size());
			computeGnMuN(gnmun,gn);

			size_t counter = 0;
			size_t n = size_t((params.omega2 - params.omega1)/params.deltaOmega); 
			if (result.size()==0) result.resize(n);
			RealType offset = params_.Eg;
			std::cerr<<"gn[0]="<<gn[0]<<" gn[5]="<<gn[5]<<"\n";
			for (RealType omega = params.omega1;omega <params.omega2;omega+=params.deltaOmega) {
				RealType x = (omega+offset-params_.b)*params_.oneOverA;
				
				RealType den = (x>1.0 || x<-1.0) ? 0.0 : sqrt(1.0 - x*x);
				RealType tmp = (fabs(den)>1e-6) ? calcF(x,gnmun)/den : 0.0;
				std::pair<RealType,RealType> p(omega,tmp);
				result[counter++] = p;

				if (counter>=result.size()) break;
				//std::cout<<omega<<" "<<real(res)<<" "<<imag(res)<<"\n";
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
			throw std::runtime_error("iOfOmega: unimplemented\n");
		}

	private:

		RealType calcF(const RealType& x,const std::vector<RealType>& gnmn) const
		{
			RealType sum = 0.5*gnmn[0];
			for (size_t i=1;i<gnmn.size();i++) sum += gnmn[i]*chebyshev_(i,x);
			return 2.0*sum;
		}
		
		void computeGnMuN( std::vector<RealType>& gnmn,std::vector<RealType>& gn) const
		{
			for (size_t i=0;i<gnmn.size();i++) gnmn[i] = moments_[i] * gn[i];
		}
		
		void initKernel(std::vector<RealType>& gn,
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

		void initKernelJackson(std::vector<RealType>& gn) const
		{
			size_t nPlus1 = gn.size()+1;
			RealType cot1 = 1.0/tan(M_PI/nPlus1);
			for (size_t i=0;i<gn.size();i++) {
				gn[i] = (nPlus1-i)*cos(M_PI*i/nPlus1)+sin(M_PI*i/nPlus1)*cot1;
				gn[i] /= nPlus1;
			}
		}

		void initKernelLorentz(std::vector<RealType>& gn,
		                       const RealType& lambda) const
		{
			RealType nreal = gn.size();
			RealType sinhlambda = sinh(lambda);
			for (size_t i=0;i<gn.size();i++) {
				gn[i] = sinh(lambda*(1-i/nreal))/sinhlambda;
			}
		}
		ProgressIndicator progress_;
		std::vector<RealType> moments_;
		ParametersType params_;
		ChebyshevFunction<RealType> chebyshev_;
	}; // class ChebyshevSerializer
	
	template<typename RealType,typename VectorType>
	const std::string ChebyshevSerializer<RealType,VectorType>::stringMarker_ =
	                                               "#ChebyshevSerializerMarker";
} // namespace PsimagLite 
/*@}*/
#endif  //CHEBYSHEV_SERIALIZER_H
