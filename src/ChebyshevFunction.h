
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

/*! \file ChebyshevFunction.h
 *
 * The Chebyshev function
 *
 */

#ifndef CHEBYSHEV_FUNCTION_H
#define CHEBYSHEV_FUNCTION_H
#include <iostream>
#include "TypeToString.h"

namespace PsimagLite {
template<typename RealType>
class ChebyshevFunction  {

public:

	RealType operator()(int m,const RealType& x) const
	{
		if (m==0) return 1;

		if (m==1) return x;

		if (m&1) {
			int p=(m-1)/2;
			return (2*this->operator()(p,x)*this->operator()(p+1,x)-x);
		}

		int pp = m/2;
		RealType tmp=this->operator()(pp,x);
		return (2*tmp*tmp-1);
	}
}; // class ChebyshevFunction

} // namespace PsimagLite 
/*@}*/
#endif  //CHEBYSHEV_FUNCTION_H

