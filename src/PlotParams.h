
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
 * Parameters for an I(omega) plot
 * This is a structure, don't add members functions (except ctor)
 *
 */
namespace PsimagLite {
template<typename RealType>
struct PlotParams {
	PlotParams(const RealType& wbegin,const RealType& wend,const RealType& wstep,const RealType& wdelta)
	: omega1(wbegin),omega2(wend),deltaOmega(wstep),delta(wdelta) { }
	RealType omega1;
	RealType omega2;
	RealType deltaOmega;
	RealType delta;
};

#ifndef PLOT_PARAMS_H
#define PLOT_PARAMS_H

} // namespace PsimagLite 
/*@}*/
#endif  //PLOT_PARAMS_H
