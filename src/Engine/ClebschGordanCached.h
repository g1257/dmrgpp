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

/** \ingroup DMRG */
/*@{*/

/*! \file ClebschGordanCached.h
 *
 *
 *
 */
#ifndef CLEBSCH_GORDANCACHED_H
#define CLEBSCH_GORDANCACHED_H

#include "ClebschGordan.h"
#include "ProgressIndicator.h"

namespace Dmrg {

template<typename FieldType>
class ClebschGordanCached {
	typedef ClebschGordan<FieldType> ClebschGordanType;
	typedef typename ClebschGordanType::PairType PairType;

public:

	ClebschGordanCached(SizeType jmax)
	    :  UNDEFINED_VALUE(-1000),
	      jmax_(jmax),
	      max2_(((jmax_-1)*(jmax_+2))/2+1),max22_(max2_*max2_),
	      data_(max22_*jmax_*2,UNDEFINED_VALUE),cgObject_(2)
	{
		copies_ = 0;
		init(jmax,2);
	}

	void init(SizeType jmax,SizeType nfactorials)
	{
		jmax_=jmax;
		max2_=((jmax_-1)*(jmax_+2))/2+1;
		max22_=max2_*max2_;
		data_.resize(max22_*jmax_*2,UNDEFINED_VALUE);
		cgObject_.init(nfactorials);

		copies_++;
		if (copies_>2) {
			std::cerr<<"WARNING: ClebschGordanCached has ";
			std::cerr<<copies_<<" copies.\n";
		}
	}

	FieldType operator()(const PairType& jm,const PairType& jm1,const PairType& jm2)
	{
		if (!checkCg(jm,jm1,jm2)) return 0;

		SizeType index1 = calcSubIndex(jm1);
		SizeType index2 = calcSubIndex(jm2);
		SizeType jmin=0;
		if (jm1.first>jm2.first) jmin = jm1.first-jm2.first;
		else jmin = jm2.first-jm1.first;
		SizeType x = calcIndex(index1,index2)+(jm.first-jmin)*max22_;
		if (data_[x]== UNDEFINED_VALUE) {
			data_[x]=cgObject_(jm,jm1,jm2);
		}
		return data_[x];
	}

private:
	SizeType calcSubIndex(const PairType& jm) const
	{
		if (jm.second==0) return jm.first;
		SizeType x = 2*jmax_-1-jm.second;
		x *= jm.second;
		x /= 2;
		x += jm.first;
		if (x>=max2_)
			throw PsimagLite::RuntimeError("problem calcSubIndex\n");
		return x;
	}

	SizeType calcIndex(SizeType i1,SizeType i2) const
	{
		if (i1>=max2_ || i2>=max2_) throw PsimagLite::RuntimeError("problem\n");
		return i1+i2*max2_;
	}

	bool checkCg(const PairType& jm,const PairType& jm1,const PairType& jm2) const
	{
		if (!checkCg1(jm.first,jm1.first,jm2.first)) return false;
		int m=calcM(jm.first,jm1,jm2);
		if (m<0) return false;
		if (SizeType(m)!=jm.second) return false;
		return true;

	}

	bool checkCg1(SizeType j,SizeType j1,SizeType j2) const
	{
		if (j>j1+j2) return false;
		SizeType jmin=0;
		if (j1<j2) jmin = j2-j1;
		else jmin  = j1-j2;
		if (j<jmin) return false;

		return true;
	}

	int calcM(SizeType j,const PairType& jm1,const PairType& jm2) const
	{
		int x = jm1.first+jm2.first-j;
		if (x%2!=0) return -1;
		x =x/2;
		if (x<0 || jm1.second+jm2.second<SizeType(x)) return -1;
		return jm1.second+jm2.second-x;
	}

	static SizeType copies_;
	int UNDEFINED_VALUE;
	SizeType jmax_;
	SizeType max2_,max22_;
	typename PsimagLite::Vector<FieldType>::Type data_;
	ClebschGordanType cgObject_;
}; // class ClebschGordanCached

template<typename FieldType>
SizeType ClebschGordanCached<FieldType>::copies_=0;

} // namespace Dmrg
/*@}*/
#endif

