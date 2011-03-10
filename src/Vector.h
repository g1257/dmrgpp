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

// temporary hack until Vector is used:
namespace std {
	template<class X>
	inline X operator*(std::vector<X> const &v,std::vector<X> const &w)
	{
		X result=0;
		for (size_t i=0;i<v.size();i++) result += v[i]*conj(w[i]);
		return result;
	}

	template<class X>
	std::ostream &operator<<(std::ostream &s,std::vector<X> const &v)
	{
		s<<v.size()<<"\n";
		for (size_t i=0;i<v.size();i++) s<<i<<" "<<v[i]<<"\n";
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
}// namespace PsimagLite
#endif // PSIVECTOR_H_

