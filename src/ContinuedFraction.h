
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
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
/** \ingroup DMRG */
/*@{*/

/*! \file ContinuedFraction.h
 *
 * A continued fraction as explained in, e.g.,
 * E. Dagotto, Rev. Mod. Phys., 66, 763, (2004).
 */

#ifndef CONTINUED_FRACTION_H
#define CONTINUED_FRACTION_H
#include <iostream>
#include "Complex.h"
#include "TypeToString.h"
#include "ProgressIndicator.h"
#include "Random48.h"

namespace PsimagLite {
	template<
		typename RealType,
    	typename TridiagonalMatrixType_>
	class ContinuedFraction  {
	public:
		typedef TridiagonalMatrixType_ TridiagonalMatrixType;
		typedef typename std::complex<RealType> ComplexType;
		typedef typename TridiagonalMatrixType::value_type FieldType;
		typedef Matrix<FieldType> MatrixType;
		typedef std::vector<std::pair<RealType,ComplexType> > PlotDataType;

		ContinuedFraction(
				const TridiagonalMatrixType& ab,
				const RealType& Eg,
				RealType weight = 1)
			: progress_("ContinuedFraction",0),ab_(ab),Eg_(Eg),weight_(weight)
		{
			if (weight_==0) return;
			MatrixType T;
			ab_.buildDenseMatrix(T);
			eigs_.resize(T.n_row());
			diag(T,eigs_,'V');
			intensity_.resize(T.n_row());
			for (size_t i=0;i<T.n_row();i++) {
				intensity_[i]= T(i,0)*T(i,0);
			}
		}

		template<typename IoInputType>
		ContinuedFraction(IoInputType& io)
		: progress_("ContinuedFraction",0),ab_(io)
		{
			io.readline(weight_,"#CFWeight=");
			io.readline(Eg_,"#CFEnergy=");
			io.read(eigs_,"#CFEigs");
			io.read(intensity_,"#CFIntensities");
		}
		
		template<typename IoOutputType>
		void save(IoOutputType& io) const
		{
			ab_.save(io);

			std::string s = "#CFWeight=" + typeToString(weight_);
			io.printline(s);

			s = "#CFEnergy=" + typeToString(Eg_);
			io.printline(s);

			io.printVector(eigs_,"#CFEigs");
			io.printVector(intensity_,"#CFIntensities");
		}

		void plot(
				PlotDataType& result,
				const RealType& omega1,
				const RealType& omega2,
				const RealType& deltaOmega,
				const RealType& delta) const
		{
			size_t counter = 0;
			size_t n = (omega2 - omega1)/deltaOmega + 1; 
			if (result.size()==0) result.resize(n);
			for (RealType omega = omega1;omega <omega2;omega+=deltaOmega) {
				ComplexType z(omega,delta);
				ComplexType res = iOfOmega(z,Eg_);
				std::pair<RealType,ComplexType> p(omega,res);
				result[counter++] = p;
				//std::cout<<omega<<" "<<real(res)<<" "<<imag(res)<<"\n";
			}
		} 

		ComplexType iOfOmega(const ComplexType& z,RealType offset) const

		{
			if (weight_==0) return ComplexType(0,0);

			ComplexType sum = 0;
			for (size_t l=0;l<intensity_.size();l++)
				sum +=intensity_[l]/(z-eigs_[l]+offset);

			return sum*weight_;
		}

	private:
		ProgressIndicator progress_;
		const TridiagonalMatrixType& ab_;
		RealType Eg_;
		RealType weight_;
		std::vector<RealType> eigs_;
		std::vector<RealType> intensity_;
	}; // class ContinuedFraction
} // namespace PsimagLite 
/*@}*/
#endif  //CONTINUED_FRACTION_H
