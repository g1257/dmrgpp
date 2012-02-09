
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

/*! \file ChebyshevFunctionExplicit.h
 *
 * The Chebyshev function with caching
 * 
 */

#ifndef CHEBYSHEV_F_EXPLICIT_H
#define CHEBYSHEV_F_EXPLICIT_H
#include "ChebyshevFunction.h"


namespace PsimagLite {
	template<typename RealType>
	class ChebyshevFunctionExplicit  {

	public:
		
		RealType operator()(int m,const RealType& x) const
		{
			return cos(m*acos(x));
		}
	}; // class ChebyshevFunctionExplicit
} // namespace PsimagLite 
/*@}*/
#endif  //CHEBYSHEV_F_EXPLICIT_H
