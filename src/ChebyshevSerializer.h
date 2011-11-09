
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
#include "Complex.h"
#include "TypeToString.h"
#include "ProgressIndicator.h"
#include "Random48.h"
#include <stdexcept>
#include "ParametersForSolver.h"

// #include "PlotParams.h"

namespace PsimagLite {
	template<
		typename RealType,
    	typename VectorType_>
	class ChebyshevSerializer  {
	public:

		typedef VectorType_ VectorType;
		typedef typename std::complex<RealType> ComplexType;
		typedef typename VectorType::value_type FieldType;
		typedef Matrix<FieldType> MatrixType;
		typedef std::vector<std::pair<RealType,ComplexType> > PlotDataType;
		typedef PlotParams<RealType> PlotParamsType;
		typedef ParametersForSolver<RealType> ParametersType;

		//! ab.a contains the even moments
		//! ab.b contains the odd moments
		//! weight contains < phi|A^\dagger A |phi>
		//! isign can be omitted
		ChebyshevSerializer(const VectorType& ab,const ParametersType& params)
		: progress_("ChebyshevSerializer",0),
		  moments_(ab),
		  Eg_(params.Eg),
		  oneOverA_((2-1e-3)/(params.eMax-params.eMin)),
		  b_((params.eMax+params.eMin)/2)
		{
			//gatherEvenAndOdd(ab);
		}

// 		ChebyshevSerializer() 
// 		: progress_("ChebyshevSerializer",0), ab_(),Eg_(0),weight_(0),isign_(1)
// 		{ }

		template<typename IoInputType>
		ChebyshevSerializer(IoInputType& io)
		: progress_("ChebyshevSerializer",0)
		{
 			io.readline(Eg_,"#ChebyshevEnergy=");
			io.readline(oneOverA_,"#ChebyshevOneOverA=");
			io.readline(b_,"#ChebyshevB");
			io.read(moments_,"#ChebyshevMoments");
		}

		template<typename IoOutputType>
		void save(IoOutputType& io) const
		{
			std::string s("#ChebyshevEnergy=");
			s += ttos(Eg_);
			io.printline(s);
			s = "#ChebyshevOneOverA=" + ttos(oneOverA_);
			io.printline(s);
			s = "#ChebyshevB=" + ttos(b_);
			io.printline(s);

			io.printVector(moments_,"#ChebyshevMoments");
		}

// 		void set(const VectorType& ab,
// 		         RealType Eg,
// 		         RealType weight = 0,
// 		         int isign = 0)
// 		{
// 			moments_=ab;
// 			//gatherEvenAndOdd(ab);
// 		}

		void plot(PlotDataType& result,const PlotParamsType& params) const
		{
			std::vector<RealType> gn(moments_.size());
			initKernelJackson(gn,moments_.size());

			std::vector<RealType> gnmun(gn.size());
			computeGnMuN(gnmun,gn);

			size_t counter = 0;
			size_t n = size_t((params.omega2 - params.omega1)/params.deltaOmega); 
			if (result.size()==0) result.resize(n);
			for (RealType omega = params.omega1;omega <params.omega2;omega+=params.deltaOmega) {
				RealType x = (omega-b_)*oneOverA_;
				RealType den = sqrt(1.0 - x*x);
				std::pair<RealType,RealType> p(omega,calcF(x,gnmun)/den);
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
		ComplexType iOfOmega(const ComplexType& z,RealType offset,int isign) const

		{
			throw std::runtime_error("iOfOmega: unimplemented\n");
		}

	private:
		
// 		void gatherEvenAndOdd(const VectorType& ab)
// 		{
// 			// need to get momentsEven, momentsOdd, a, and b
// 			size_t n = ab.size();
// 			std::vector<RealType> moments(2*n);
// 			for (size_t i=0;i<n;i++) {
// 				moments[2*i] = ab.a(i);
// 				moments[2*i+1] = ab.b(i);
// 			}
// 		}
		
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

		void initKernelJackson(std::vector<RealType>& gn) const
		{
			size_t nPlus1 = gn.size()+1;
			RealType cot1 = 1.0/tan(M_PI/nPlus1);
			for (size_t i=0;i<gn.size();i++) {
				gn[i] = (nPlus1-i)*cos(M_PI*i/nPlus1)+sin(M_PI*i/nPlus1)*cot1;
				gn[i] /= nPlus1;
			}
		}

		ProgressIndicator progress_;
		std::vector<RealType> moments_;
		RealType Eg_;
		RealType oneOverA_,b_;
	}; // class ChebyshevSerializer
} // namespace PsimagLite 
/*@}*/
#endif  //CHEBYSHEV_SERIALIZER_H
