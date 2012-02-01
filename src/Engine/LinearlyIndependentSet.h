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

/*! \file LinearlyIndependentSet
 *
 *
 */
#ifndef LINEARLY_INDEPENDENT_SET
#define LINEARLY_INDEPENDENT_SET

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"
//#include "Vector.h"
#include "SparseVector.h"

namespace Dmrg {

// FIXME: MOVE ELSEWHERE:
template<typename RealType>
bool isAlmostZero(const RealType& x,double eps = 1e-16)
{
	return (fabs(x)<eps);
}

// FIXME: MOVE ELSEWHERE:
template<typename RealType>
bool isAlmostZero(const std::complex<RealType>& x,double eps = 1e-16)
{
	return (fabs(real(x)*real(x)+imag(x)*imag(x))<eps);
}

template<typename RealType,typename SparseMatrixType>
class LinearlyIndependentSet {

	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef SparseVector<ComplexOrRealType> SparseVectorType;

public:

	LinearlyIndependentSet(size_t rank)
	: transform_(rank,rank),counter_(0),row_(0)
	{}

	~LinearlyIndependentSet()
	{
		deallocate();
	}

	void pushNew(SparseVectorType& v2)
	{
		v2.sort();
		RealType norma = PsimagLite::norm(v2);
		if (isAlmostZero(norma,1e-8)) return;

		v2 *= (1.0/norma);

		for (size_t i=0;i<e_.size();i++) {
			const SparseVectorType& v3 = *e_[i];
			if (v3 == v2) return;
			v2 *= (-1);
			if (v3 == v2) return;
		}
		e_.push_back(&v2);
		fill(v2);
	}

	void push(SparseVectorType& v2)
	{
		v2.sort();
		RealType norma = PsimagLite::norm(v2);
		if (isAlmostZero(norma,1e-8)) return;

		v2 *= (1.0/norma);

		if (e_.size()==0) {
//			vecs_.push_back(&v2);
			SparseVectorType* u = new SparseVectorType(v2);
			e_.push_back(u);
			fill(v2);
			return;
		}

		SparseVectorType* u = new SparseVectorType(v2);
		ComplexOrRealType x = (v2*v2);
		for (size_t i=0;i<e_.size();i++) {
			SparseVectorType tmp = ((*e_[i])*v2)*(*e_[i]);
			x -= (tmp*v2);
			(*u) -= tmp;
		}

		if (std::norm(x)<1e-6) return;

		u->sort();
		norma = PsimagLite::norm(*u);
		(*u) *= (1.0/norma);

//		vecs_.push_back(&v2);
		e_.push_back(u);
		//printVector(v2,"ADDING: ");
		fill(v2);
	}

	void clear()
	{
		deallocate();
	}

	size_t size() const { return e_.size(); }

	const SparseMatrixType& transform()
	{
		return transform_;
	}

private:


	void fill(SparseVectorType& vref)
	{
		std::cerr<<__FILE__<<" push vecs.size="<<e_.size()<<"\n";
		//std::cerr<<vref<<"\n";
		size_t i = row_;
		assert(row_<transform_.rank());
		transform_.setRow(i,counter_);
		for (size_t k=0;k<vref.indices();k++) {
			size_t j = vref.index(k);
			ComplexOrRealType val = vref.value(k); //vecs_[i][j];
			if (isAlmostZero(val,1e-8)) continue;
			transform_.pushCol(j);
			transform_.pushValue(val);
			counter_++;
		}
		row_++;
		transform_.setRow(transform_.rank(),counter_);
	}

	void deallocate()
	{
		for (size_t i=0;i<e_.size();i++)
			delete e_[i];
		e_.clear();
	}

	void printVector(const SparseVectorType& v,const std::string& label) const
	{
		std::cout<<label;
		for (size_t i=0;i<v.size();i++)
			std::cout<<v[i]<<" ";
		std::cout<<"\n";
	}

	SparseMatrixType transform_;
	size_t counter_;
	size_t row_;
//	std::vector<SparseVectorType*> vecs_;
	std::vector<SparseVectorType*> e_;

}; // class LinearlyIndependentSet

//template<typename T1,typename T2>
//std::ostream& operator<<(std::ostream& os,const LinearlyIndependentSet<T1,T2>& ri)
//{
//	ri.print(os);
//	return os;
//}

} // namespace Dmrg 

/*@}*/
#endif // LINEARLY_INDEPENDENT_SET
