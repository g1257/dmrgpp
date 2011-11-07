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

namespace std {

	template<class T1,class T2>
	ostream &operator<<(std::ostream &os,const pair<T1,T2>& p)
	{
		os<<p.first<<" "<<p.second;
		return os;
	}
} // namespace std

//! Utility functions that are still needed
namespace utils {
	bool neighbors(size_t i1,size_t i2)
	{
		return (i1-i2==1 || i2-i1==1);
	}
	
	template<typename SomeType>
	void reorder(std::vector<SomeType> &v,std::vector<size_t> const &permutation)
	{
		std::vector<SomeType> tmpVector(v.size());
		for (size_t i=0;i<v.size();i++) tmpVector[i]=v[permutation[i]]; 
		v = tmpVector;
	}
	
	template<typename SomeType>
	void reorder(PsimagLite::Matrix<SomeType>& v,std::vector<size_t> const &permutation)
	{
		PsimagLite::Matrix<SomeType> tmpVector(v.n_row(),v.n_col());
		for (size_t i=0;i<v.n_row();i++) 
			for (size_t j=0;j<v.n_col();j++)
				tmpVector(i,j)=v(permutation[i],permutation[j]); 
		v = tmpVector;
	}
	
	//! A = B union C
        template<typename Block>
        void blockUnion(Block &A,Block const &B,Block const &C)
        {
		A=B;
		for (size_t i=0;i<C.size();i++) A.push_back(C[i]);
        }
        
        template<typename SomeType>
	void truncateVector(std::vector<SomeType> &v,std::vector<SomeType> const &removedIndices)
	{
		std::vector<SomeType> tmpVector;
		for (size_t i=0;i<v.size();i++) {
			if (PsimagLite::isInVector(removedIndices,i)>=0) continue;
			tmpVector.push_back(v[i]);
		}
		v=tmpVector;
	}
        
        template<class T>
	void truncate(PsimagLite::Matrix<T> &A,std::vector<size_t> const &removed,bool rowOption)
	{
		size_t j;
		int x=removed.size();
		if (x<=0) return;
		size_t nrow = A.n_row();
		size_t ncol = A.n_col();
		
		size_t n = ncol;
		if (rowOption)  n = nrow;
		
		if (int(n)<=x) {
			std::cerr<<"psimag::truncate: n="<<n<<" must be larger than x="<<x<<" rowoption="<<rowOption<<"\n";
			throw std::runtime_error("psimag::truncated\n");
		}
		
		std::vector<int> remap(n);
		
		
		//! find remapping
		j=0;
		for (size_t i=0;i<n;i++) {
			remap[i] = -1;
			if (PsimagLite::isInVector(removed,i)>=0) continue;
			remap[i]=j;
			j++;
		}
		if (j!=n-x) throw std::runtime_error("truncate: PsimagLite::Matrix is throwing...\n");
		
		//! truncate
		if (rowOption) {
			PsimagLite::Matrix<T> B(nrow-x,ncol);
			for (size_t i=0;i<ncol;i++) {
				for (j=0;j<nrow;j++) {
					if (remap[j]<0) continue;
					B(remap[j],i)=A(j,i);
				}
			}
			A=B; 
		} else {
			PsimagLite::Matrix<T> B(nrow,ncol-x);
			for (size_t i=0;i<nrow;i++) {
				for (j=0;j<ncol;j++) {
					if (remap[j]<0) continue;
					B(i,remap[j])=A(i,j);
				}
			}
			A=B; 
		}
	}

	
} //namespace utils
/*@}*/
#endif


