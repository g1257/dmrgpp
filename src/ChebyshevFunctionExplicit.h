
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
			
			switch (m) {
				case 0:
					return 1;
				case 1:
					return x;
			}
			
			RealType x1 = x*x;
			switch (m) {
				case 2:
					return 2*x1-1;
				case 3:
					return 4*x1*x-3*x;
				case 4:
					return 8*x1*x1-8*x1+1;
				case 5:
					return 16*x1*x1*x1-20*x1*x+5*x;
			}
			
			RealType x2 = x1*x1;
			switch (m) {
				case 6:
					return 32*x2*x1-48*x2+18*x1-1;
				case 7:
					return 64*x2*x1*x-112*x2*x+56*x1*x-7*x;
			}
			
			if (m&1) {
				int p=(m-1)/2;
				return (2*this->operator()(p,x)*this->operator()(p+1,x)-x);
			}

			int pp = m/2;
			RealType tmp=this->operator()(pp,x);
			return (2*tmp*tmp-1);
		}
	}; // class ChebyshevFunctionExplicit
} // namespace PsimagLite 
/*@}*/
#endif  //CHEBYSHEV_F_EXPLICIT_H
