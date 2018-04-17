/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file JmPairs.h
 *
 *  This is a "vector" of (2j,m+j) pairs.
 *  Repeated entries are stored only once
 *  Provides a transparent access as if it were a normal vector
 *
 */
#ifndef JMPAIRS_HEADER_H
#define JMPAIRS_HEADER_H

#include <algorithm>
#include "Utils.h"
#include "PsimagLite.h"

namespace Dmrg {

template<typename PairType_>
class JmPairs {

	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	typedef PairType_ PairType;
	typedef PairType value_type;

	PairType operator[](SizeType alpha) const
	{
		return jmPairs_[indices_[alpha]];
	}

	JmPairs<PairType>& operator=(const typename PsimagLite::Vector<PairType>::Type& jmpairs)
	{
		jmPairs_.clear();
		indices_.clear();
		for (SizeType i=0;i<jmpairs.size();i++) {
			int x = PsimagLite::isInVector(jmPairs_,jmpairs[i]);
			if (x<0) {
				jmPairs_.push_back(jmpairs[i]);
				x=jmPairs_.size()-1;
			}
			indices_.push_back(x);
		}
		return *this;
	}

	//! indices_[alpha] = jm
	void push(const PairType& jm,SizeType)
	{
		int x = PsimagLite::isInVector(jmPairs_,jm);

		if (x<0) {
			jmPairs_.push_back(jm);
			x=jmPairs_.size()-1;
		}

		indices_.push_back(x);
	}

	void clear()
	{
		jmPairs_.clear();
		indices_.clear();
	}

	void reorder(const VectorSizeType& permutation)
	{
		utils::reorder(indices_,permutation);
	}

	void truncate(const VectorSizeType& removedIndices)
	{
		utils::truncateVector(indices_,removedIndices);
		VectorSizeType unusedPairs;
		findUnusedJmPairs(unusedPairs);
		removeUnusedPairs(unusedPairs);
	}

	template<typename Op>
	void maxFirst(SizeType& maxvalue)
	{
		Op f;
		for (SizeType i=0;i<jmPairs_.size();i++) {
			if (f(jmPairs_[i].first,maxvalue)) {
				maxvalue=jmPairs_[i].first;
			}
		}
	}

	SizeType size() const { return indices_.size(); }

	void resize(SizeType) { } // does nothing, safely

	template<typename IoOutputter>
	void write(IoOutputter& io,
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsOutputLike<IoOutputter>::True, int>::Type = 0) const
	{
		io.write(jmPairs_,"#su2JmPairs");
		io.write(indices_,"#su2JmIndices");
	}

	template<typename IoInputter>
	void read(IoInputter& io,
	          typename PsimagLite::EnableIf<
	          PsimagLite::IsInputLike<IoInputter>::True, int>::Type = 0)
	{
		io.read(jmPairs_,"#su2JmPairs");
		io.read(indices_,"#su2JmIndices");
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                JmPairs<PairType> jmPairs)
	{
		for (SizeType i=0;i<jmPairs.size();i++)
			os<<"jmPair["<<i<<"]="<<jmPairs[i]<<"\n";
		return os;
	}

private:

	void findUnusedJmPairs(VectorSizeType& unusedPairs)
	{
		for (SizeType i=0;i<jmPairs_.size();i++)
			if (isUnusedPair(i)) unusedPairs.push_back(i);
	}

	void removeUnusedPairs(const VectorSizeType& unusedPairs)
	{
		SizeType counter=0;
		VectorSizeType neworder(jmPairs_.size());
		typename PsimagLite::Vector<PairType>::Type tmpVector(jmPairs_.size() -
		                                                      unusedPairs.size());

		for (SizeType i=0;i<jmPairs_.size();i++) {
			if (PsimagLite::isInVector(unusedPairs,i)>=0) continue;
			tmpVector[counter]=jmPairs_[i];
			neworder[i]=counter;
			counter++;
		}

		jmPairs_=tmpVector;
		VectorSizeType tmpVector2(indices_.size());
		for (SizeType i=0;i<indices_.size();i++)
			tmpVector2[i]=neworder[indices_[i]];
		indices_=tmpVector2;
	}

	bool isUnusedPair(SizeType ind)
	{
		for (SizeType i=0;i<indices_.size();i++)
			if (indices_[i]==ind) return false;
		return true;
	}

	typename PsimagLite::Vector<PairType>::Type jmPairs_;
	VectorSizeType indices_;
}; // JmPairs
} // namespace Dmrg
/*@}*/
#endif

