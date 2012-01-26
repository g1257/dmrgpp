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

/*! \file ReflectionItem
 *
 *
 */
#ifndef REFLECTION_ITEM_H
#define REFLECTION_ITEM_H

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"

namespace Dmrg {

// FIXME: MOVE ELSEWHERE:
template<typename RealType>
bool isAlmostZero(const RealType& x)
{
	return (fabs(x)<1e-6);
}

// FIXME: MOVE ELSEWHERE:
template<typename RealType>
bool isAlmostZero(const std::complex<RealType>& x)
{
	return (fabs(real(x)*real(x)+imag(x)*imag(x))<1e-6);
}

template<typename RealType,typename ComplexOrRealType>
class ReflectionItem {

public:

	enum { DIAGONAL,PLUS,MINUS};

	ReflectionItem(size_t i)
	: type_(DIAGONAL),i_(i),vec_(0)
	{}

	ReflectionItem(size_t type1,size_t i,const std::vector<ComplexOrRealType>& v)
	: type_(type1),i_(i),vec_(v)
	{}

	size_t type() const { return type_; }

	void postFix()
	{
		if (type_==DIAGONAL) return;
		setVector();
		int x = getUniqueIndex();
		if (x>=0 && size_t(x)==i_) type_ = DIAGONAL;
	}

	template<typename SomeSparseMatrix>
	void setTransformPlus(SomeSparseMatrix& transform,size_t& counter) const
	{
		assert(type_!=MINUS);
		switch(type_) {
		case DIAGONAL:
			transform.pushCol(i_);
			transform.pushValue(1.0);
			counter++;
			break;
		case PLUS:
			for (size_t k=0;k<vec_.size();k++) {
				ComplexOrRealType val =vec_[k];
				if (isAlmostZero(val)) continue;
				transform.pushCol(k);
				transform.pushValue(val);
				counter++;
			}
			break;
		}
	}

	template<typename SomeSparseMatrix>
	void setTransformMinus(SomeSparseMatrix& transform,size_t& counter) const
	{
		assert(type_==MINUS);
		for (size_t k=0;k<vec_.size();k++) {
			ComplexOrRealType val = vec_[k];
			if (isAlmostZero(val)) continue;
			transform.pushCol(k);
			transform.pushValue(val);
			counter++;
		}
	}

	bool operator==(const ReflectionItem<RealType,ComplexOrRealType>& item2) const
	{
		if (type_>item2.type_) return (item2==*this);

		if (type_==DIAGONAL) {
			if (item2.type_==DIAGONAL) {
				return (i_==item2.i_);
			}
			return false;
		}

		std::vector<ComplexOrRealType> v3=(-1.0)*vec_;
		return (equalV(vec_,item2.vec_) || equalV(v3,item2.vec_));
	}

	void print(std::ostream& os) const
	{
		std::string s = typeToString();
		s += " i=" + ttos(i_);
		os<<s<<" ";
		if (type_==DIAGONAL) {
			os<<"\n";
			return;
		}
		for (size_t i=0;i<vec_.size();i++) {
			os<<vec_[i]<<" ";
		}
		os<<"\n";

	}

private:

	std::string typeToString() const
	{
		switch(type_) {
		case DIAGONAL:
			return "DIAGONAL";
		case PLUS:
			return "PLUS";
		case MINUS:
			return "MINUS";
		}
		assert(false);
		return "UNDEFINED";
	}

	bool equalV(const std::vector<ComplexOrRealType>& v1,const std::vector<ComplexOrRealType>& v2) const
	{
		for (size_t i=0;i<v1.size();i++) {
			ComplexOrRealType x = v1[i]-v2[i];
			if (!isAlmostZero(x)) return false;
		}
		return true;
	}

	int getUniqueIndex() const
	{
		assert(type_!=DIAGONAL);
		size_t j=0;
		bool seenBefore = false;
		for (size_t i=0;i<vec_.size();i++) {
			if (isAlmostZero(vec_[i])) continue;
			if (seenBefore) return -1;
			seenBefore=true;
			j=i;
		}
		if (!seenBefore) return -1;
		return j;
	}

	void setVector()
	{
		RealType plusOrMinus = (type_==PLUS) ? 1 : -1;
		vec_[i_]+=plusOrMinus;
		RealType norm = PsimagLite::norm(vec_);
		assert(norm>1e-6);
		vec_/=norm;
	}

	size_t type_,i_;
	std::vector<ComplexOrRealType> vec_;
}; // class ReflectionItem

template<typename T1,typename T2>
std::ostream& operator<<(std::ostream& os,const ReflectionItem<T1,T2>& ri)
{
	ri.print(os);
	return os;
}

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_ITEM_H
