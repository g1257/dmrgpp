/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

#ifndef PSICOMPLEX_H_
#define PSICOMPLEX_H_

#include <complex>
#include "../loki/TypeTraits.h"

namespace PsimagLite {

inline double real(double t) { return t; }

inline double real(const std::complex<double>& t)
{
	return std::real(t);
}

inline double imag(double) { return 0.0; }

inline double imag(const std::complex<double>& t)
{
	return std::imag(t);
}

inline double conj(double t) { return t; }

template<typename T>
inline std::complex<T> conj(const std::complex<T>& t) { return std::conj(t); }

inline double norm(double t)
{
	return fabs(t);
}

template<typename T>
inline double norm(const std::complex<T>& t)
{
	return std::norm(t);
}

}

namespace std {
inline std::complex<double> operator*(int x,const std::complex<double>& y)
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

template<typename T>
struct IsNumber {
	enum {True = (IsComplexNumber<T>::True || Loki::TypeTraits<T>::isArith)};
};
} // namespace PsimagLite

#endif // PSICOMPLEX_H_

