// BEGIN LICENSE BLOCK
/*
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
// END LICENSE BLOCK
#ifndef PSICOMPLEX_H_
#define PSICOMPLEX_H_

#include <complex>
namespace std {
	double real(double t) { return t; }

	double imag(double t) { return 0.0; }

	double conj(double t) { return t; }
	
	double norm(double t)
	{
		return fabs(t);
	}

	std::complex<double> operator*(int x,const std::complex<double>& y)
	{
		return std::complex<double>(real(y)*x,imag(y)*x);
	}
} // namespace std

namespace PsimagLite {

	template<typename ComplexOrRealType>
	class Real {
	public:
		typedef ComplexOrRealType Type;
	};

	template<typename RealType>
	class Real<std::complex<RealType> > {
	public:
		typedef RealType Type;
	};

	template<typename T>
	class IsComplexNumber {
	public:
		enum { True = false};
	};

	template<typename T>
	class IsComplexNumber<std::complex<T> > {
	public:
		enum { True = true};
	};
} // namespace PsimagLite

#endif // PSICOMPLEX_H_
