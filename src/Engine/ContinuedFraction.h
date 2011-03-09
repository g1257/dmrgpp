
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef CONTINUED_FRACTION_H 
#define CONTINUED_FRACTION_H
#include <iostream>
#include "ProgressIndicator.h"
#include "Random48.h"
#include "ProgramGlobals.h"

namespace Dmrg {
	template<
		typename RealType,
    	typename TridiagonalMatrixType>
	class ContinuedFraction  {
	public:
		
		typedef typename std::complex<RealType> ComplexType;
		typedef typename TridiagonalMatrixType::FieldType FieldType;
		typedef PsimagLite::Matrix<FieldType> MatrixType;
		
		ContinuedFraction(
				const TridiagonalMatrixType& ab,
				const RealType& Eg,
				RealType weight = 1)
			: progress_("ContinuedFraction",0),ab_(ab),weight_(1),Eg_(Eg)
		{
			MatrixType T;
			ab_.buildDenseMatrix(T);
			intensity_.resize(T.n_row());
			for (size_t i=0;i<T.n_row();i++) {
				intensity_[i]= T(i,0)*T(i,0);
			}
			eigs_.resize(T.n_row());
			diag(T,eigs_,'V');
		}

		template<typename IoInputType>
		ContinuedFraction(IoInputType& io)
		: progress_("ContinuedFraction",0),ab_(io)
		{
			io.readline("#CFWeight",weight_);
			io.readline("#CFEnergy",Eg_);
		}
		
		template<typename IoOutputType>
		void save(IoOutputType& io) const
		{
			ab_.save(io);

			std::string s = "#CFWeight=" + utils::ttos(weight_);
			io.printline(s);

			std::string s = "#CFEnergy=" + utils::ttos(Eg_);
			io.printline(s);
		}

		void plot(
				const RealType& omega1,
				const RealType& omega2,
				const RealType& deltaOmega,
				const RealType& delta) const
		{
			for (RealType omega = omega1;omega <omega2;omega+=deltaOmega) {
				ComplexType z(omega,delta);
				ComplexType res = calcIntensity(z);
				std::cout<<omega<<" "<<real(res)<<" "<<imag(res)<<"\n";
			}
		} 

		ComplexType iOfOmega(const ComplexType& z,RealType offset = Eg_) const

		{
			ComplexType sum = 0;
			for (size_t l=0;l<S.n_row();l++)
				sum +=intensity_[i]/(z-eigs[l]+offset);

			return sum*weight_;
		}
		
	private:
		PsimagLite::ProgressIndicator progress_;
		const TridiagonalMatrixType& ab_;
		RealType Eg_;
		RealType weight_;
		std::vector<RealType> eigs_;
		std::vector<RealType> intensity_;
	}; // class ContinuedFraction
} // namespace Dmrg

#endif 
