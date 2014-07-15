
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
#ifndef PLOT_PARAMS_H
#define PLOT_PARAMS_H

namespace PsimagLite {
template<typename RealType>
struct PlotParams {

	PlotParams(RealType wbegin,
	           RealType wend,
	           RealType wstep,
	           RealType wdelta,
	           RealType beta1,
	           SizeType numberOfMatsubaras1)
	: omega1(wbegin),omega2(wend),deltaOmega(wstep),delta(wdelta),
	  beta(beta1),numberOfMatsubaras(numberOfMatsubaras1)
	{}

	RealType omega1;
	RealType omega2;
	RealType deltaOmega;
	RealType delta;
	RealType beta;
	SizeType numberOfMatsubaras;
};
} // namespace PsimagLite
/*@}*/
#endif  //PLOT_PARAMS_H

