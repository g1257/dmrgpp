/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."

*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************

*/
#ifndef UTILS_HEADER_H
#define UTILS_HEADER_H

#include "Vector.h"
#include "CrsMatrix.h"

namespace std {

template<class T1,class T2>
ostream &operator<<(std::ostream &os,const pair<T1,T2>& p)
{
	os<<p.first<<" "<<p.second;
	return os;
}
} // namespace std

// Utility functions that are still needed
namespace utils {

struct UnixPathSeparator
{
	bool operator()(char ch) const
	{
		return ch == '/';
	}
};

PsimagLite::String basename(PsimagLite::String pathname);

PsimagLite::String pathPrepend(PsimagLite::String,PsimagLite::String);

template<template<typename,typename> class SomeVectorTemplate,
         typename SomeAllocator1Type,
         typename SomeAllocator2Type,
         typename T>
typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<
SomeVectorTemplate<T,SomeAllocator1Type> >::True,void>::Type
reorder(SomeVectorTemplate<T,SomeAllocator1Type>& v,
        const SomeVectorTemplate<SizeType,SomeAllocator2Type>& permutation)
{
	SomeVectorTemplate<T,SomeAllocator1Type> tmpVector(v.size());
	for (SizeType i=0;i<v.size();i++) tmpVector[i]=v[permutation[i]];
	v = tmpVector;
}

template<typename SomeType>
void reorder(PsimagLite::Matrix<SomeType>& v,
             const PsimagLite::Vector<SizeType>::Type& permutation)
{
	PsimagLite::Matrix<SomeType> tmpVector(v.n_row(),v.n_col());
	for (SizeType i=0;i<v.n_row();i++)
		for (SizeType j=0;j<v.n_col();j++)
			tmpVector(i,j)=v(permutation[i],permutation[j]);
	v = tmpVector;
}

//! A = B union C
template<typename Block>
void blockUnion(Block &A,Block const &B,Block const &C)
{
	A=B;
	for (SizeType i=0;i<C.size();i++) A.push_back(C[i]);
}

template<typename SomeVectorType>
typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,void>::Type
truncateVector(SomeVectorType& v,
               const SomeVectorType& removedIndices)
{
	SomeVectorType tmpVector;
	for (SizeType i=0;i<v.size();i++) {
		if (PsimagLite::isInVector(removedIndices,i)>=0) continue;
		tmpVector.push_back(v[i]);
	}
	v=tmpVector;
}

template<typename SomeVectorType>
static
typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,void>::Type
fillFermionicSigns(SomeVectorType& fermionicSigns,
                   const typename PsimagLite::Vector<SizeType>::Type& electrons,
                   int f)
{
	typedef typename SomeVectorType::value_type ValueType;
	fermionicSigns.resize(electrons.size());
	for (SizeType i=0;i<fermionicSigns.size();i++)
		fermionicSigns[i]= (electrons[i]%2==0) ? 1.0 : static_cast<ValueType>(f);
}

SizeType exactDivision(SizeType,SizeType);

SizeType bitSizeOfInteger(SizeType);

SizeType powUint(SizeType, SizeType);

} //namespace utils
/*@}*/
#endif

