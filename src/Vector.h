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
#ifndef PSIVECTOR_H_
#define PSIVECTOR_H_

#include <vector>
#include <iostream>
#include <stdexcept>
#include "Complex.h"

// temporary hack until Vector is used:
namespace std {
template<class X>
inline X operator*(const std::vector<X>& v,const std::vector<X>& w)
{
	X result=0;
	for (size_t i=0;i<v.size();i++) result += v[i]*conj(w[i]);
	return result;
}

template<typename T1,typename T2>
inline std::vector<T2> operator*(const std::vector<std::vector<T1> >& v1,const std::vector<T2>& v2)
{
	std::vector<T2> v3(v2.size());
	for (size_t i=0;i<v3.size();i++) {
		v3[i] = 0;
		for (size_t j=0;i<v2.size();j++)
			v3[i] += v1[i][j] * v2[j];
	}
	return v3;
}

template<typename T1,typename T2>
inline std::vector<T2> operator*(const T1& v1,const std::vector<T2>& v2)
{
	std::vector<T2> v3(v2.size());
	for (size_t i=0;i<v3.size();i++) v3[i] = v1 * v2[i];
	return v3;
}

template<typename T1,typename T2>
inline std::vector<T2> operator*(const std::vector<T2>& v2,const T1& v1)
{
	return v1*v2;
}

template<typename T>
std::vector<T> conj(std::vector<T>& v)
{
	std::vector<T> w(v.size());
	for (size_t i=0;i<v.size();i++) w[i]=std::conj(v[i]);
	return w;
}

template<typename T>
T scalarProduct(const std::vector<T>& v1, const std::vector<T>& v2)
{
	T result = 0.0;
	for(size_t i=0; i < v2.size(); i++)
		result += std::conj(v1[i]) * v2[i];
	return result;
}

template<typename FieldType>
inline std::vector<FieldType> operator+=(std::vector<FieldType>& v,const std::vector<FieldType>& w)
{
	for (size_t i=0;i<w.size();i++) v[i] += w[i];
	return v;
}

template<typename FieldType>
inline std::vector<FieldType> operator-=(std::vector<FieldType>& v,const std::vector<FieldType>& w)
{
	for (size_t i=0;i<w.size();i++) v[i] -= w[i];
	return v;
}

template<typename T1,typename T2>
inline std::vector<T1> operator*=(std::vector<T1>& v,const T2& t2)
{
	for (size_t i=0;i<v.size();i++) v[i] *= t2;
	return v;
}

template<typename T1,typename T2>
inline std::vector<T1> operator/=(std::vector<T1>& v,const T2& t2)
{
	for (size_t i=0;i<v.size();i++) v[i] /= t2;
	return v;
}

template<typename T1,typename T2>
inline std::vector<T1> operator+(const std::vector<T1>& v1,const std::vector<T2>& v2)
{
	std::vector<T1> v3(v1.size());
	for (size_t i=0;i<v1.size();i++) v3[i] = v1[i] + v2[i];
	return v3;
}

template<typename T1,typename T2>
inline std::vector<T1> operator-(const std::vector<T1>& v1,const std::vector<T2>& v2)
{
	std::vector<T1> v3(v1.size());
	for (size_t i=0;i<v1.size();i++) v3[i] = v1[i] - v2[i];
	return v3;
}

template<class X>
std::ostream &operator<<(std::ostream &s,const std::vector<X>& v)
{
	s<<v.size()<<"\n";
	for (size_t i=0;i<v.size();i++) s<<v[i]<<"\n";
	return s;
}

template<typename FieldType>
inline std::istream& operator>>(std::istream& is,std::vector<FieldType>& v)
{
	int xsize = 0;
	is>>xsize;
	if (xsize<0) throw std::runtime_error(">> vector: size is negative\n");
	v.resize(xsize);
	for (size_t i=0;i<size_t(xsize);i++) {
		is>>v[i];
	}
	return is;
}
} // namespace std 

namespace PsimagLite {
	// FIXME: write proper class here
	template<typename T>
	class  Vector  {
	public:
		typedef std::vector<T> Type;
	}; // class Vector

	// change this when using PsimagLite::Vector: 
	template<class T>
	void vectorPrint(const std::vector<T>& v,char const *name,std::ostream &s)
	{
		for (size_t i=0;i<v.size();i++) s<<name<<"["<<i<<"]="<<v[i]<<std::endl;
	} 
	
	template<class X>
	X norm(const std::vector<X>& v)
	{
		return sqrt(v*v);
	}

	template<class X>
	X norm(const std::vector<std::complex<X> >& v)
	{
		std::complex<X> x = v*v;
		if (fabs(imag(x))>1e-5) throw std::runtime_error("Norm isn't real\n");
		return sqrt(real(x));
	}

	template<typename X,typename RandomType>
	void randomizeVector(typename Vector<typename RandomType::value_type>::Type& v,const X& a,const X& b,const RandomType& r)
	{
		for (size_t i=0;i<v.size();i++) v[i] = a + b*r.random();
	}
	
	template<typename X,typename Y>
	int isInVector(const std::vector<X>& natBasis,Y const &v)
	{
		typename std::vector<X>::const_iterator x = find(natBasis.begin(),natBasis.end(),v);
		if (x==natBasis.end()) return -1;
		return x-natBasis.begin();
		
	}

	template<typename T>
	void split(std::vector<T>& v,const char* s1,char sep)
	{
		std::string buffer = "";
		std::string s(s1);
		T tmp;
		for (size_t i=0;i<s.length();i++) {
			if (s[i]==sep) {
				if (buffer=="") continue;
				std::istringstream buffer2(buffer);
				buffer2 >> tmp;
				buffer="";
				v.push_back(tmp);
			} else {
				buffer += s[i];
			}
		}
		if (buffer=="") return;
		std::istringstream buffer2(buffer);
		buffer2 >> tmp;
		buffer="";
		v.push_back(tmp);
	}
}// namespace PsimagLite
#endif // PSIVECTOR_H_

