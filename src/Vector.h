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
	inline X operator*(std::vector<X> const &v,std::vector<X> const &w)
	{
		X result=0;
		for (size_t i=0;i<v.size();i++) result += v[i]*conj(w[i]);
		return result;
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

	template<class X>
	std::ostream &operator<<(std::ostream &s,std::vector<X> const &v)
	{
		s<<v.size()<<"\n";
		for (size_t i=0;i<v.size();i++) s<<v[i]<<"\n";
		return s;
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
	class  Vector : public std::vector<T> {
	public:
		typedef T ValueType;
	}; // class Vector

	// change this when using PsimagLite::Vector: 
	template<class T>
	void vectorPrint(std::vector<T> const &v,char const *name,std::ostream &s)
	{
		for (size_t i=0;i<v.size();i++) s<<name<<"["<<i<<"]="<<v[i]<<std::endl;
	} 
	
	template<class X>
	X norm(std::vector<X> const &v)
	{
		return sqrt(v*v);
	}

	template<class X>
	X norm(std::vector<std::complex<X> > const &v)
	{
		std::complex<X> x = v*v;
		if (fabs(imag(x))>1e-5) throw std::runtime_error("Norm isn't real\n");
		return sqrt(real(x));
	}

	template<typename X,typename RandomType>
	void randomizeVector(std::vector<typename RandomType::value_type>& v,const X& a,const X& b,const RandomType& r)
	{
		for (size_t i=0;i<v.size();i++) v[i] = a + b*r.random();
	}
	
	template<typename X,typename Y>
	int isInVector(std::vector<X> const &natBasis,Y const &v)
	{
// 		if (natBasis.size()==0) return -1;
// 		for (size_t ii=0;ii<natBasis.size();ii++) if (natBasis[ii]==v) return ii;
// 		return -1;
		typename std::vector<X>::const_iterator x = find(natBasis.begin(),natBasis.end(),v);
		if (x==natBasis.end()) return -1;
		return x-natBasis.begin();
		
	}
}// namespace PsimagLite
#endif // PSIVECTOR_H_

