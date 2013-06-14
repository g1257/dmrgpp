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
/** \ingroup PsimagLite */
/*@{*/

/*! \file Vector.h
 *
 *
 */

#ifndef PSIVECTOR_H_
#define PSIVECTOR_H_

#include <vector>
#include <iostream>
#include <stdexcept>
#include "Complex.h"
#include "AllocatorCpu.h"
#include "String.h"

// temporary hack until Vector is used:
namespace std {

template<typename T1,typename T2>
inline std::istream& operator>>(std::istream& is,std::pair<T1,T2>& v)
{
	is>>v.first;
	is>>v.second;
	return is;
}

template<typename X,typename A>
inline X operator*(const std::vector<X,A>& v,const std::vector<X,A>& w)
{
	X result=0;
	for (SizeType i=0;i<v.size();i++) result += v[i]*conj(w[i]);
	return result;
}

template<typename T1,typename T2,typename A,typename AA>
inline std::vector<T2,A> operator*(const std::vector<std::vector<T1,A>,AA>& v1,
                                 const std::vector<T2,A>& v2)
{
	std::vector<T2,A> v3(v2.size());
	for (SizeType i=0;i<v3.size();i++) {
		v3[i] = 0;
		for (SizeType j=0;i<v2.size();j++)
			v3[i] += v1[i][j] * v2[j];
	}
	return v3;
}

template<typename T1,typename T2,typename A>
inline std::vector<T2,A> operator*(const T1& v1,const std::vector<T2,A>& v2)
{
	std::vector<T2,A> v3(v2.size());
	for (SizeType i=0;i<v3.size();i++) v3[i] = v1 * v2[i];
	return v3;
}

template<typename T1,typename T2,typename A>
inline std::vector<T2,A> operator*(const std::vector<T2,A>& v2,const T1& v1)
{
	return v1*v2;
}

template<typename T,typename A>
std::vector<T,A> conj(std::vector<T,A>& v)
{
	std::vector<T,A> w(v.size());
	for (SizeType i=0;i<v.size();i++) w[i]=std::conj(v[i]);
	return w;
}

template<typename T,typename A>
T scalarProduct(const std::vector<T,A>& v1, const std::vector<T,A>& v2)
{
	T result = 0.0;
	for(SizeType i=0; i < v2.size(); i++)
		result += std::conj(v1[i]) * v2[i];
	return result;
}

template<typename FieldType,typename A>
inline std::vector<FieldType,A> operator+=(std::vector<FieldType,A>& v,const std::vector<FieldType,A>& w)
{
	for (SizeType i=0;i<w.size();i++) v[i] += w[i];
	return v;
}

template<typename FieldType,typename A>
inline std::vector<FieldType,A> operator-=(std::vector<FieldType,A>& v,const std::vector<FieldType,A>& w)
{
	for (SizeType i=0;i<w.size();i++) v[i] -= w[i];
	return v;
}

template<typename T1,typename T2,typename A>
inline std::vector<T1,A> operator*=(std::vector<T1,A>& v,const T2& t2)
{
	for (SizeType i=0;i<v.size();i++) v[i] *= t2;
	return v;
}

template<typename T1,typename T2,typename A>
inline std::vector<T1,A> operator/=(std::vector<T1,A>& v,const T2& t2)
{
	for (SizeType i=0;i<v.size();i++) v[i] /= t2;
	return v;
}

template<typename T1,typename T2,typename A1,typename A2>
inline std::vector<T1,A1> operator+(const std::vector<T1,A1>& v1,const std::vector<T2,A2>& v2)
{
	std::vector<T1,A1> v3(v1.size());
	for (SizeType i=0;i<v1.size();i++) v3[i] = v1[i] + v2[i];
	return v3;
}

template<typename T1,typename T2,typename A>
inline std::vector<T1,A> operator-(const std::vector<T1,A>& v1,const std::vector<T2,A>& v2)
{
	std::vector<T1,A> v3(v1.size());
	for (SizeType i=0;i<v1.size();i++) v3[i] = v1[i] - v2[i];
	return v3;
}

template<class X,typename A>
std::ostream &operator<<(std::ostream &s,const std::vector<X,A>& v)
{
	s<<v.size()<<"\n";
	for (SizeType i=0;i<v.size();i++) s<<v[i]<<"\n";
	return s;
}

template<typename FieldType,typename A>
inline std::istream& operator>>(std::istream& is,std::vector<FieldType,A>& v)
{
	int xsize = 0;
	is>>xsize;
	if (xsize<0)
		throw PsimagLite::RuntimeError(">> vector: size is negative\n");
	v.resize(xsize);
	for (SizeType i=0;i<SizeType(xsize);i++) {
		is>>v[i];
	}
	return is;
}
} // namespace std 

namespace PsimagLite {

	template<typename T>
	class  Vector  {
	public:
		typedef std::vector<T,typename Allocator<T>::Type> Type;
	}; // class Vector

	template<>
	class  Vector<bool>  {
	public:
		typedef std::vector<bool,Allocator<bool>::Type> Type;
	}; // class Vector

	template<typename T>
	class IsVectorLike {
	public:
		enum {True = false};
	};

	template<typename T>
	class IsVectorLike<std::vector<T,typename Allocator<T>::Type> > {
	public:
	        enum {True = true};
	};

	// change this when using PsimagLite::Vector: 
	template<class T,typename A>
	void vectorPrint(const std::vector<T,A>& v,char const *name,std::ostream &s)
	{
		for (SizeType i=0;i<v.size();i++) s<<name<<"["<<i<<"]="<<v[i]<<std::endl;
	} 
	
	template<class X,typename A>
	X norm(const std::vector<X,A>& v)
	{
		return sqrt(v*v);
	}

	template<class X,typename A>
	X norm(const std::vector<std::complex<X>,A>& v)
	{
		std::complex<X> x = v*v;
		if (fabs(imag(x))>1e-5) throw RuntimeError("Norm isn't real\n");
		return sqrt(real(x));
	}

	template<typename X,typename RandomType>
	void randomizeVector(typename Vector<typename RandomType::value_type>::Type& v,
	                     const X& a,
	                     const X& b,
	                     const RandomType& r)
	{
		for (SizeType i=0;i<v.size();i++) v[i] = a + b*r.random();
	}
	
	template<typename X,typename Y,typename A>
	int isInVector(const std::vector<X,A>& natBasis,Y const &v)
	{
		typename std::vector<X,A>::const_iterator x = find(natBasis.begin(),natBasis.end(),v);
		if (x==natBasis.end()) return -1;
		return x-natBasis.begin();
		
	}

	template<typename SomeVectorType>
	typename EnableIf<IsVectorLike<SomeVectorType>::True,
	                  typename SomeVectorType::value_type>::Type
	sum(SomeVectorType& v)
	{
		typename SomeVectorType::value_type tmp = 0;
		for (size_t i=0;i<v.size();i++) {
			tmp += v[i];
		}
		return tmp;
	}

	template<typename T,typename A>
	void split(std::vector<T,A>& v,const char* s1,char sep)
	{
		String buffer = "";
		String s(s1);
		T tmp;
		for (SizeType i=0;i<s.length();i++) {
			if (s[i]==sep) {
				if (buffer=="") continue;
				PsimagLite::IstringStream buffer2(buffer);
				buffer2 >> tmp;
				buffer="";
				v.push_back(tmp);
			} else {
				buffer += s[i];
			}
		}
		if (buffer=="") return;
		PsimagLite::IstringStream buffer2(buffer);
		buffer2 >> tmp;
		buffer="";
		v.push_back(tmp);
	}
}// namespace PsimagLite

/*@}*/
#endif // PSIVECTOR_H_

