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
#include "Vector.h"

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

template<typename RealType,typename VectorType>
class LinearlyIndependentSet {

	typedef typename VectorType::value_type ComplexOrRealType;

public:

	LinearlyIndependentSet(size_t rank)
		: rank_(rank)
	{}

	void push(const VectorType& v)
	{
		RealType norma = PsimagLite::norm(v);
		if (isAlmostZero(norma)) return;

		VectorType v2;
		v2 = (1.0/norma)*v;

		if (vecs_.size()==0) {
			vecs_.push_back(v2);
			printVector(v2,"ADDING: ");
			e_.push_back(v2);
			return;
		}

		VectorType u = v2;
		for (size_t i=0;i<vecs_.size();i++) {
			u -= (e_[i]*v2)*e_[i];
		}

		ComplexOrRealType x = (u*v2);

		if (std::norm(x)<1e-6) return;

		norma = PsimagLite::norm(u);
		u /= norma;

		vecs_.push_back(v2);
		e_.push_back(u);
		printVector(v2,"ADDING: ");
	}

	size_t size() const { return vecs_.size(); }

	template<typename SomeSparseMatrixType>
	void fill(SomeSparseMatrixType& s)
	{
		s.resize(vecs_.size());
		size_t counter=0;
		for (size_t i=0;i<vecs_.size();i++) {
			s.setRow(i,counter);
			for (size_t j=0;j<vecs_[i].size();j++) {
				ComplexOrRealType val = vecs_[i][j];
				if (isAlmostZero(val)) continue;
				s.pushCol(j);
				s.pushValue(val);
				counter++;
			}
		}
		s.setRow(s.rank(),counter);
	}

private:

	void printVector(const VectorType& v,const std::string& label) const
	{
		std::cout<<label;
		for (size_t i=0;i<v.size();i++)
			std::cout<<v[i]<<" ";
		std::cout<<"\n";
	}

	size_t rank_;
	std::vector<VectorType> vecs_;
	std::vector<VectorType> e_;

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
