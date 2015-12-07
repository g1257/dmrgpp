// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
// END LICENSE BLOCK
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

PsimagLite::String basename(PsimagLite::String pathname)
{
	return PsimagLite::String(std::find_if(pathname.rbegin(),
	                                       pathname.rend(),
	                                       UnixPathSeparator()).base(), pathname.end());
}

template<template<typename,typename> class SomeVectorTemplate,
         typename SomeAllocator1Type,
         typename SomeAllocator2Type,
         typename T>
typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorTemplate<T,SomeAllocator1Type> >::True,void>::Type
reorder(SomeVectorTemplate<T,SomeAllocator1Type>& v,
        const SomeVectorTemplate<SizeType,SomeAllocator2Type>& permutation)
{
	SomeVectorTemplate<T,SomeAllocator1Type> tmpVector(v.size());
	for (SizeType i=0;i<v.size();i++) tmpVector[i]=v[permutation[i]];
	v = tmpVector;
}

template<typename SomeType>
void reorder(PsimagLite::Matrix<SomeType>& v,const PsimagLite::Vector<SizeType>::Type& permutation)
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

template<class T>
void truncate(PsimagLite::Matrix<T> &A,const PsimagLite::Vector<SizeType>::Type& removed,bool rowOption)
{
	SizeType j;
	int x=removed.size();
	if (x<=0) return;
	SizeType nrow = A.n_row();
	SizeType ncol = A.n_col();

	SizeType n = ncol;
	if (rowOption)  n = nrow;

	if (int(n)<=x) {
		std::cerr<<"psimag::truncate: n="<<n<<" must be larger than x="<<x<<" rowoption="<<rowOption<<"\n";
		throw PsimagLite::RuntimeError("psimag::truncated\n");
	}

	PsimagLite::Vector<int>::Type remap(n);


	//! find remapping
	j=0;
	for (SizeType i=0;i<n;i++) {
		remap[i] = -1;
		if (PsimagLite::isInVector(removed,i)>=0) continue;
		remap[i]=j;
		j++;
	}
	if (j!=n-x) throw PsimagLite::RuntimeError("truncate: PsimagLite::Matrix is throwing...\n");

	//! truncate
	if (rowOption) {
		PsimagLite::Matrix<T> B(nrow-x,ncol);
		for (SizeType i=0;i<ncol;i++) {
			for (j=0;j<nrow;j++) {
				if (remap[j]<0) continue;
				B(remap[j],i)=A(j,i);
			}
		}
		A=B;
	} else {
		PsimagLite::Matrix<T> B(nrow,ncol-x);
		for (SizeType i=0;i<nrow;i++) {
			for (j=0;j<ncol;j++) {
				if (remap[j]<0) continue;
				B(i,remap[j])=A(i,j);
			}
		}
		A=B;
	}
}

template<class T>
void truncate(PsimagLite::CrsMatrix<T> &A,const PsimagLite::Vector<SizeType>::Type& removed,bool rowOption)
{
	if (rowOption) { // unimplemented
		assert(false);
		throw PsimagLite::RuntimeError("truncate: rowoption must not be set\n");
	}

	SizeType x=removed.size();
	if (x==0) return;

	SizeType nrow = A.row();

	SizeType n = nrow;

	assert(n>x);

	PsimagLite::Vector<int>::Type remap(n);

	//! find remapping
	SizeType j=0;
	for (SizeType i=0;i<n;i++) {
		remap[i] = -1;
		if (PsimagLite::isInVector(removed,i)>=0) continue;
		remap[i]=j;
		j++;
	}
	assert(j==n-x);

	//! truncate
	PsimagLite::CrsMatrix<T> B(nrow,nrow-x);
	SizeType counter = 0;
	for (SizeType i=0;i<nrow;i++) {
		B.setRow(i,counter);
		for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
			j = A.getCol(k);
			if (remap[j]<0) continue;
			B.pushCol(remap[j]);
			B.pushValue(A.getValue(k));
			counter++;
			//B(i,remap[j])=A(i,j);
		}
	}
	B.setRow(nrow,counter);
	B.checkValidity();
	A=B;

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

SizeType exactDivision(SizeType a,SizeType b)
{
	SizeType c = static_cast<SizeType>(a/b);
	if (c * b != a)
		throw PsimagLite::RuntimeError("exactDivision expected\n");

	return c;
}

SizeType log2OfInteger(SizeType x)
{
	SizeType counter = 0;
	while (x) {
		counter++;
		x >>= 1;
	}

	return counter;
}

} //namespace utils
/*@}*/
#endif


