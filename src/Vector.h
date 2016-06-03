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
#include "../loki/TypeTraits.h"

namespace PsimagLite {

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

} // namespace PsimagLite

namespace std {

template<typename T1,typename T2>
istream& operator>>(istream& is,pair<T1,T2>& v)
{
	is>>v.first;
	is>>v.second;
	return is;
}

template<typename X,typename A>
X operator*(const vector<X,A>& v,const vector<X,A>& w)
{
	X result=0;
	for (SizeType i=0;i<v.size();i++) result += v[i]*conj(w[i]);
	return result;
}

template<typename T1,typename T2,typename A,typename AA>
vector<T2,A> operator*(const vector<vector<T1,A>,AA>& v1,
                       const vector<T2,A>& v2)
{
	vector<T2,A> v3(v2.size());
	for (SizeType i=0;i<v3.size();i++) {
		v3[i] = 0;
		for (SizeType j=0;i<v2.size();j++)
			v3[i] += v1[i][j] * v2[j];
	}
	return v3;
}

// closure

struct ClosureOperations {

	enum {OP_PLUS,OP_MINUS,OP_MULT,OP_DIVIDE,OP_CONJ};
};

template<typename T1, typename T2,int type>
class ClosureOperator {

public:

	ClosureOperator(const T1& v1,const T2& v2)
	    : r1(v1),r2(v2)
	{}

	const T1& r1;
	const T2& r2;
};

template<typename T>
struct IsClosureLike {

	enum {True = false};
};

template<typename T1, typename T2,int type>
struct IsClosureLike<ClosureOperator<T1,T2,type> > {

	enum {True = true};
};

// vector * scalar
template<typename T1,typename T2,typename A>
ClosureOperator<T1,vector<T2,A>,ClosureOperations::OP_MULT> operator*(const T1& v1,
                                                                      const vector<T2,A>& v2)
{
	return ClosureOperator<T1,vector<T2,A>,ClosureOperations::OP_MULT >(v1,v2);
}

// vector * scalar
template<typename T1,typename T2>
typename PsimagLite::EnableIf<IsClosureLike<T1>::True || IsClosureLike<T2>::True,
ClosureOperator<T1,T2,ClosureOperations::OP_MULT> >::Type
operator*(const T1& v1,const T2& v2)
{
	return ClosureOperator<T1,T2,ClosureOperations::OP_MULT>(v1,v2);
}

template<typename T1,typename T2,typename A>
typename PsimagLite::EnableIf<PsimagLite::IsNumber<T1>::True,void>::Type
operator<=(vector<T2,A>& v,
                const ClosureOperator<T1,vector<T2,A>,ClosureOperations::OP_MULT>& c)
{
	v = c.r2;
	for (SizeType i=0;i<v.size();i++) v[i] *= c.r1;
}

template<typename T1,typename T2,typename A>
typename PsimagLite::EnableIf<PsimagLite::IsNumber<T1>::True && PsimagLite::IsNumber<T2>::True,
void>::Type operator<=(vector<T2,A>& v,
                const ClosureOperator<T1,
                ClosureOperator<T2,std::vector<T2,A>,ClosureOperations::OP_MULT>,
                ClosureOperations::OP_MULT>& c)
{
	v = c.r2.r2;
	T2 tmp = c.r1*c.r2.r1;
	for (SizeType i=0;i<v.size();i++) v[i] *= tmp;
}

template<typename T1,typename T2,typename A>
ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT> operator*(const vector<T2,A>& v2,
                                                                       const T1& v1)
{
	return v1*v2;
}

// vector + vector
template<typename T1,typename T2,typename A1, typename A2>
ClosureOperator<vector<T1,A1>,vector<T2,A2>,ClosureOperations::OP_PLUS>
operator+(const vector<T1,A1>& v1,const vector<T2,A2>& v2)
{
	return ClosureOperator<vector<T1,A1>,vector<T2,A2>,ClosureOperations::OP_PLUS>(v1,v2);
}

template<typename T1,typename T2>
typename PsimagLite::EnableIf<IsClosureLike<T1>::True || IsClosureLike<T2>::True,
ClosureOperator<T1,T2,ClosureOperations::OP_PLUS> >::Type
operator+(const T1& v1,const T2& v2)
{
	return ClosureOperator<T1,T2,ClosureOperations::OP_PLUS>(v1,v2);
}

template<typename T1,typename T2,typename A>
void operator<=(vector<T1,A>& v,
                const ClosureOperator<vector<T1,A>,vector<T2,A>,ClosureOperations::OP_PLUS>& c)
{
	v.resize(c.r1.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1[i] + c.r2[i];
}

template<typename T1,typename T2,typename A>
void operator<=(vector<T2,A>& v,
                const ClosureOperator<vector<T2,A>,
                ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_PLUS>& c)
{
	v.resize(c.r1.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1[i] + c.r2.r1*c.r2.r2[i];
}

template<typename T1,typename T2,typename A1, typename A2>
void operator<=(vector<T1,A1>& v,
                const ClosureOperator<
                ClosureOperator<T1,vector<T2,A2>, ClosureOperations::OP_MULT>,
                vector<T1,A1>,
                ClosureOperations::OP_PLUS>& c)
{
	v.resize(c.r2.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1.r1*c.r1.r2[i] + c.r2[i];
}

template<typename T1,typename T2,typename A>
void operator<=(vector<T2,A>& v,
                const ClosureOperator<T1,
                ClosureOperator<vector<T2,A>,vector<T2,A>, ClosureOperations::OP_PLUS>,
                ClosureOperations::OP_MULT>& c)
{
	v.resize(c.r2.r2.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1*(c.r2.r1[i] + c.r2.r2[i]);
}

template<typename T1,typename T2,typename A>
void operator<=(vector<T2,A>& v,
                const ClosureOperator<
                ClosureOperator<
                ClosureOperator<ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_PLUS>,
                ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_PLUS>,
                ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_PLUS>& c)
{
	v.resize(c.r2.r2.size());
	T1 m1 = c.r2.r1;
	T1 m2 = c.r1.r2.r1;
	T1 m3 = c.r1.r1.r2.r1;
	T1 m4 = c.r1.r1.r1.r1;
	const vector<T2,A>& k1 = c.r2.r2;
	const vector<T2,A>& k2 = c.r1.r2.r2;
	const vector<T2,A>& k3 = c.r1.r1.r2.r2;
	const vector<T2,A>& k4 = c.r1.r1.r1.r2;
	for (SizeType i=0;i<v.size();i++)
		v[i] = m1*k1[i] + m2*k2[i] + m3*k3[i] + m4*k4[i];
}

// vector - vector

template<typename T1,typename T2,typename A1, typename A2>
ClosureOperator<vector<T1,A1>,vector<T2,A2>,ClosureOperations::OP_MINUS>
operator-(const vector<T1,A1>& v1,const vector<T2,A2>& v2)
{
	return ClosureOperator<vector<T1,A1>,vector<T2,A2>,ClosureOperations::OP_MINUS>(v1,v2);
}

template<typename T1,typename T2>
typename PsimagLite::EnableIf<IsClosureLike<T1>::True || IsClosureLike<T2>::True,
ClosureOperator<T1,T2,ClosureOperations::OP_MINUS> >::Type
operator-(const T1& v1,const T2& v2)
{
	return ClosureOperator<T1,T2,ClosureOperations::OP_MINUS>(v1,v2);
}

template<typename T,typename A>
void operator<=(vector<T,A>& v,
                const ClosureOperator<vector<T,A>,vector<T,A>,ClosureOperations::OP_MINUS>& c)
{
	v.resize(c.r1.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1[i] - c.r2[i];
}

template<typename T1,typename T2,typename A1, typename A2>
void operator<=(vector<T1,A1>& v,
                const ClosureOperator<vector<T1,A1>,
                ClosureOperator<T1,vector<T2,A2>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_MINUS>& c)
{
	v.resize(c.r1.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1[i] - c.r2.r1*c.r2.r2[i];
}

template<typename T1,typename T2,typename A>
void operator<=(vector<T2,A>& v,
                const ClosureOperator<T1,
                ClosureOperator<vector<T2,A>,vector<T2,A>, ClosureOperations::OP_MINUS>,
                ClosureOperations::OP_MULT>& c)
{
	v.resize(c.r2.r2.size());
	for (SizeType i=0;i<v.size();i++) v[i] = c.r1*(c.r2.r1[i] - c.r2.r2[i]);
}

template<typename T1, typename T2, typename A>
void operator<=(vector<T2,A>& v,
                const ClosureOperator<
                ClosureOperator<std::vector<T2,A>,
                ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_MINUS>,
                ClosureOperator<T1,vector<T2,A>, ClosureOperations::OP_MULT>,
                ClosureOperations::OP_PLUS>& c)
{
	T1 m2 = c.r1.r2.r1;
	T1 m3 = c.r2.r1;
	const vector<T2,A>& k1 = c.r1.r1;
	const vector<T2,A>& k2 = c.r1.r2.r2;
	const vector<T2,A>& k3 = c.r2.r2;

	v.resize(k1.size());
	for (SizeType i=0;i<v.size();i++)
		v[i] = k1[i] - m2*k2[i] + m3*k3[i];
}

// operator+=
template<typename FieldType,typename A>
vector<FieldType,A> operator+=(vector<FieldType,A>& v,
                               const vector<FieldType,A>& w)
{
	for (SizeType i=0;i<w.size();i++) v[i] += w[i];
	return v;
}

template<typename T1,typename T2,typename A>
vector<T2,A> operator+=(vector<T2,A>& v,
                       const ClosureOperator<T1,vector<T2,A>,ClosureOperations::OP_MULT>& w)
{
	for (SizeType i=0;i<v.size();i++) v[i] += w.r1*w.r2[i];
	return v;
}

// operator-=
template<typename FieldType,typename A>
vector<FieldType,A> operator-=(vector<FieldType,A>& v,const vector<FieldType,A>& w)
{
	for (SizeType i=0;i<w.size();i++) v[i] -= w[i];
	return v;
}

template<typename T1,typename T2,typename A>
vector<T2,A> operator-=(vector<T2,A>& v,
                       const ClosureOperator<T1,vector<T2,A>,ClosureOperations::OP_MULT>& w)
{
	for (SizeType i=0;i<v.size();i++) v[i] -= w.r1*w.r2[i];
	return v;
}

// operator*=
template<typename T1,typename T2,typename A>
vector<T1,A> operator*=(vector<T1,A>& v,
                        const T2& t2)
{
	for (SizeType i=0;i<v.size();i++) v[i] *= t2;
	return v;
}

template<typename T1,typename T2,typename A>
vector<T1,A> operator/=(vector<T1,A>& v,
                        const T2& t2)
{
	for (SizeType i=0;i<v.size();i++) v[i] /= t2;
	return v;
}

// end of closure

template<typename T,typename A>
T scalarProduct(const vector<T,A>& v1, const vector<T,A>& v2)
{
	T result = 0.0;
	for (SizeType i=0; i < v2.size(); i++)
		result += conj(v1[i]) * v2[i];
	return result;
}

template<class X,typename A>
ostream &operator<<(ostream &s,const vector<X,A>& v)
{
	s<<v.size()<<"\n";
	for (SizeType i=0;i<v.size();i++) s<<v[i]<<"\n";
	return s;
}

template<typename X,typename Y,typename A>
ostream &operator<<(ostream &s,
                    const vector<pair<X,Y>,A>& v)
{
	s<<v.size()<<"\n";
	for (SizeType i=0;i<v.size();i++) s<<v[i].first<<" "<<v[i].second<<"\n";
	return s;
}

template<typename FieldType,typename A>
istream& operator>>(istream& is,vector<FieldType,A>& v)
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

template<typename T>
class IsPairLike {
public:
	enum {True = false};
};

template<typename T1, typename T2>
class IsPairLike<std::pair<T1,T2> > {
public:
	enum {True = true};
};

}// namespace PsimagLite

/*@}*/
#endif // PSIVECTOR_H_

