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
#include "AllocatorCpu.h"

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
	enum { True = Loki::TypeTraits<T>::isArith };
};


template<typename T>
struct IsNumber {
	enum {True = (IsComplexNumber<T>::True || Loki::TypeTraits<T>::isArith)};
};

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
real(T t) { return t; }

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
real(const std::complex<T>& t)
{
	return std::real(t);
}

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
imag(T) { return 0.0; }

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
imag(const std::complex<T>& t)
{
	return std::imag(t);
}

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
conj(T t) { return t; }

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,std::complex<T> >::Type
conj(const std::complex<T>& t) { return std::conj(t); }

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
norm(T t)
{
	return fabs(t);
}

template<typename T>
typename EnableIf<Loki::TypeTraits<T>::isFloat,T>::Type
norm(const std::complex<T>& t)
{
	return std::norm(t);
}
} // namespace PsimagLite

namespace std {
template<typename T>
typename PsimagLite::EnableIf<Loki::TypeTraits<T>::isFloat,std::complex<T> >::Type
operator*(int x,const std::complex<T>& y)
{
	return std::complex<T>(real(y)*x,imag(y)*x);
}
} // namespace std

#endif // PSICOMPLEX_H_

