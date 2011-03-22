
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

/*! \file TwoContinuedFraction.h
 *
 * We need two cont. fractions for the Green's function
 * on different sites, because, you know, we need
 * c_i + c_j and also c_i - c_j
 * This class handles the composition
 */

#ifndef TWO_CONTINUED_FRACTION_H
#define TWO_CONTINUED_FRACTION_H
#include <iostream>
#include "LineMarker.h"

namespace PsimagLite {

template<typename RealType>
std::vector<std::pair<RealType,std::complex<RealType> > > operator-(
		const std::vector<std::pair<RealType,std::complex<RealType> > >& v1,
		const std::vector<std::pair<RealType,std::complex<RealType> > >& v2)
{
	std::vector<std::pair<RealType,std::complex<RealType> > > v(v1.size());
	for (size_t i=0;i<v1.size();i++) {
		v[i].first = v1[i].first;
		v[i].second = v1[i].second - v2[i].second;
	}
	return v;
}

//template<typename RealType>
//std::vector<std::pair<RealType,std::complex<RealType> > > equal1(
//		const std::vector<std::pair<RealType,std::complex<RealType> > >& v1)
//{
//	std::vector<std::pair<RealType,std::complex<RealType> > > v(v1.size());
//	for (size_t i=0;i<v1.size();i++) {
//		v[i].first = v1[i].first;
//		v[i].second = v1[i].second;
//	}
//	return v;
//}

	template<typename ContinuedFractionType>
	class TwoContinuedFraction  {
	public:
		
		typedef typename ContinuedFractionType::ComplexType ComplexType;
		typedef typename ContinuedFractionType::TridiagonalMatrixType
				TridiagonalMatrixType;
		typedef typename	TridiagonalMatrixType::value_type RealType;
		typedef typename ContinuedFractionType::MatrixType MatrixType;
		typedef typename ContinuedFractionType::PlotDataType PlotDataType;

		TwoContinuedFraction(
				const ContinuedFractionType& cf1,
				const ContinuedFractionType& cf2)
			: progress_("TwoContinuedFraction",0),
			  lmarker_("#TWOCONTINUEDFRACTION"),cf1_(cf1),cf2_(cf2)
		{
		}

		template<typename IoInputType>
		TwoContinuedFraction(IoInputType& io,size_t level = 0)
		: progress_("ContinuedFraction",0),
		  lmarker_(io,"#TWOCONTINUEDFRACTION",level),cf1_(io),cf2_(io)
		{
		}
		
		void plot(
				PlotDataType& result,
				const RealType& omega1,
				const RealType& omega2,
				const RealType& deltaOmega,
				const RealType& delta) const
		{
			PlotDataType result1;
			cf1_.plot(result1,omega1,omega2,deltaOmega,delta);

			PlotDataType result2;
			cf2_.plot(result2,omega1,omega2,deltaOmega,delta);

			result = result1 - result2;
		}
		
		template<typename IoOutputType>
		void save(IoOutputType& io) const
		{
			lmarker_.save(io);
			cf1_.save(io);
			cf2_.save(io);
		}

	private:
		ProgressIndicator progress_;
		PsimagLite::LineMarker lmarker_;
		ContinuedFractionType cf1_,cf2_;
	}; // class TwoContinuedFraction
} // namespace PsimagLite 
/*@}*/
#endif  //TWO_CONTINUED_FRACTION_H
