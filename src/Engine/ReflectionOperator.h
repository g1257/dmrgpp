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
/** \ingroup DMRG */
/*@{*/

/*! \file ReflectionOperator
 *
 *  Problems
 *  - Getting transform from stacked transform at wft
 *  - Getting bases for sys and env from stacks at checkpoint
 *  - non-existant permutations
 *
 */
#ifndef REFLECTION_OPERATOR_H
#define REFLECTION_OPERATOR_H

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"

namespace Dmrg {

class ReflectionItem {

public:

	enum { DIAGONAL,PLUS,MINUS};

	ReflectionItem(size_t ii)
		: i(ii),j(ii),type(DIAGONAL)
	{}

	ReflectionItem(size_t ii,size_t jj,size_t type1)
		: i(ii),j(jj),type(type1)
	{}

	size_t i,j,type;

}; // class ReflectionItem

bool operator==(const ReflectionItem& item1,const ReflectionItem& item2)
{
	if (item1.type!=item2.type) return false;

	if (item1.i==item2.j && item1.j==item2.i) return true;

	return (item1.i==item2.i && item1.j==item2.j);
}

template<typename LeftRightSuperType>
class ReflectionOperator {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::SparseMatrixType
			SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef ReflectionItem ItemType;
	typedef std::vector<ComplexOrRealType> VectorType;

public:

	ReflectionOperator(const LeftRightSuperType& lrs,size_t n0,bool isEnabled)
	: lrs_(lrs),
	  n0_(n0),
	  progress_("ReflectionOperator",0),
	  plusSector_(0),
	  isEnabled_(isEnabled)
	{}

	void check(const std::vector<size_t>& sectors)
	{
		if (!isEnabled_) return;
		assert(sectors.size()==1);
		size_t m = sectors[0];
		size_t offset = lrs_.super().partition(m);
		size_t total = lrs_.super().partition(m+1)-offset;
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		size_t counter = 0;
		ComplexOrRealType one = 1.0;
		s_.resize(total);
		for (size_t i=0;i<total;i++) {
			s_.setRow(i,counter);
			size_t x=0,y=0;
			pack1.unpack(x,y,lrs_.super().permutation(i+offset));

			size_t x0=0,x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			size_t y0=0,y1=0;
			pack3.unpack(y0,y1,lrs_.right().permutation(y));

			size_t xprime = pack2.pack(x0,y0,lrs_.left().permutationInverse());
			size_t yprime = pack3.pack(x1,y1,lrs_.right().permutationInverse());
			size_t iprime = pack1.pack(xprime,yprime,lrs_.super().permutationInverse());

			assert(iprime>=offset && iprime<total+offset);
			s_.pushCol(iprime-offset);
			s_.pushValue(one);
			counter++;
		}
		s_.setRow(total,counter);
		//print();
		computeTransform();
	}

	void transform(SparseMatrixType& matrixA,
		       SparseMatrixType& matrixB,
		       const SparseMatrixType& matrix) const
	 {
		assert(isEnabled_);

		 SparseMatrixType rT;
		 transposeConjugate(rT,transform_);

		 SparseMatrixType tmp;
		 multiply(tmp,matrix,rT);

		 SparseMatrixType matrix2;
		 multiply(matrix2,transform_,tmp);

		 split(matrixA,matrixB,matrix2);
	 }

	template<typename SomeVectorType>
	void setInitState(const SomeVectorType& initVector,
			  SomeVectorType& initVector1,
			  SomeVectorType& initVector2) const
	{
		assert(isEnabled_);
		size_t minusSector = initVector.size()-plusSector_;
		initVector1.resize(plusSector_);
		initVector2.resize(minusSector);
		for (size_t i=0;i<initVector.size();i++) {
			if (i<plusSector_) initVector1[i]=initVector[i];
			else  initVector2[i-plusSector_]=initVector[i];
		}

	}

	RealType setGroundState(VectorType& gs,
				const RealType& gsEnergy1,
				const VectorType& gsVector1,
				const RealType& gsEnergy2,
				const VectorType& gsVector2) const
	{
		assert(isEnabled_);
		size_t rank = gsVector1.size() + gsVector2.size();
		if (gsEnergy1<=gsEnergy2) {
			setGs(gs,gsVector1,rank,0);
			return gsEnergy1;
		}
		setGs(gs,gsVector2,rank,gsVector1.size());
		return gsEnergy2;
	}

	bool isEnabled() const { return isEnabled_; }

private:

	void setGs(VectorType& gs,const VectorType& v,size_t rank,size_t offset) const
	{
		VectorType gstmp(rank,0);

		for (size_t i=0;i<v.size();i++) {
			gstmp[i+offset]=v[i];
		}
		SparseMatrixType rT;
		transposeConjugate(rT,transform_);
		multiply(gs,rT,gstmp);
	}

	void computeTransform()
	{
		std::vector<ItemType> items;
		ItemType item(0);
		for (size_t i=0;i<s_.rank();i++) {
			for (int k=s_.getRowPtr(i);k<s_.getRowPtr(i+1);k++) {
				size_t col = s_.getCol(k);
				if (col==i) {
					item.type = ItemType::DIAGONAL;
					item.i=item.j=i;
					items.push_back(item);
					continue;
				}
				// add plus
				item.type = ItemType::PLUS;
				item.i = i;
				item.j = col;
				items.push_back(item);
				// add minus
				item.type = ItemType::MINUS;
				item.i = i;
				item.j = col;
				items.push_back(item);
			}
		}
		setTransform(items);
	}

	void setTransform(const std::vector<ItemType>& buffer2)
	{
		std::vector<ItemType> buffer;
		makeUnique(buffer,buffer2);
		transform_.resize(s_.rank());
		assert(buffer.size()==transform_.rank());
		size_t counter = 0;
		RealType oneOverSqrt2 = 1.0/sqrt(2.0);
		RealType sign = 1.0;
		size_t row = 0;
		for (size_t i=0;i<buffer.size();i++) {
			if (buffer[i].type==ItemType::MINUS) continue;
			transform_.setRow(row++,counter);
			switch(buffer[i].type) {
			case ItemType::DIAGONAL:
				transform_.pushCol(buffer[i].i);
				transform_.pushValue(1);
				counter++;
				break;
			case ItemType::PLUS:
				transform_.pushCol(buffer[i].i);
				transform_.pushValue(oneOverSqrt2);
				counter++;
				transform_.pushCol(buffer[i].j);
				transform_.pushValue(oneOverSqrt2);
				counter++;
				break;
			}
		}

		for (size_t i=0;i<buffer.size();i++) {
			if (buffer[i].type!=ItemType::MINUS) continue;
			transform_.setRow(row++,counter);
			transform_.pushCol(buffer[i].i);
			transform_.pushValue(oneOverSqrt2*sign);
			counter++;
			transform_.pushCol(buffer[i].j);
			transform_.pushValue(-oneOverSqrt2*sign);
			counter++;
		}
		transform_.setRow(transform_.rank(),counter);
	}

	void makeUnique(std::vector<ItemType>& dest,const std::vector<ItemType>& src)
	{
		size_t zeros=0;
		size_t pluses=0;
		size_t minuses=0;
		for (size_t i=0;i<src.size();i++) {
			ItemType item = src[i];
			int x =  PsimagLite::isInVector(dest,item);
			if (x>=0) continue;
//				if (item.type ==ItemType::PLUS) {
//					size_t i = item.i;
//					size_t j = item.j;
//					ItemType item2(j,i,ItemType::PLUS);
//					x = PsimagLite::isInVector(dest,item2);
//					if (x>=0) continue;
//				}
			if (item.type==ItemType::DIAGONAL) zeros++;
			if (item.type==ItemType::PLUS) pluses++;
			if (item.type==ItemType::MINUS) minuses++;

			dest.push_back(item);
		}
		std::ostringstream msg;
		msg<<pluses<<" +, "<<minuses<<" -, "<<zeros<<" zeros.";
		progress_.printline(msg,std::cout);
		plusSector_ = zeros + pluses;
	}

	void print()
	{
		PsimagLite::Matrix<ComplexOrRealType> fullm(s_.rank(),s_.rank());
		crsMatrixToFullMatrix(fullm,s_);
		std::cout<<"-------------------\n";
		std::cout<<fullm;

	}

	bool isAlmostZero(const RealType& x) const
	{
		return (fabs(x)<1e-6);
	}

	bool isAlmostZero(const std::complex<RealType>& x) const
	{
		return (fabs(real(x)*real(x)+imag(x)*imag(x))<1e-6);
	}

	void split(SparseMatrixType& matrixA,SparseMatrixType& matrixB,const SparseMatrixType& matrix) const
	{
		size_t counter = 0;
		matrixA.resize(plusSector_);
		for (size_t i=0;i<plusSector_;i++) {
			matrixA.setRow(i,counter);
			for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
				size_t col = matrix.getCol(k);
				ComplexOrRealType val = matrix.getValue(k);
				if (col<plusSector_) {
					matrixA.pushCol(col);
					matrixA.pushValue(val);
					counter++;
					continue;
				}
				if (!isAlmostZero(val)) {
					std::string s(__FILE__);
					s += " Hamiltonian has no reflection symmetry.";
					throw std::runtime_error(s.c_str());
				}
			}
		}
		matrixA.setRow(plusSector_,counter);

		size_t rank = matrix.rank();
		size_t minusSector=rank-plusSector_;
		matrixB.resize(minusSector);
		counter=0;
		for (size_t i=plusSector_;i<rank;i++) {
			matrixB.setRow(i-plusSector_,counter);
			for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
				size_t col = matrix.getCol(k);
				ComplexOrRealType val = matrix.getValue(k);
				if (col>=plusSector_) {
					matrixB.pushCol(col-plusSector_);
					matrixB.pushValue(val);
					counter++;
					continue;
				}
				if (!isAlmostZero(val)) {
					std::string s(__FILE__);
					s += " Hamiltonian has no reflection symmetry.";
					throw std::runtime_error(s.c_str());
				}
			}
		}
		matrixB.setRow(minusSector,counter);
	}

	//SparseMatrixType& operator()() const { return s_; }

	const LeftRightSuperType& lrs_;
	size_t n0_; // hilbert size of one site
	PsimagLite::ProgressIndicator progress_;
	size_t plusSector_;
	bool isEnabled_;
	SparseMatrixType s_;
	SparseMatrixType transform_;
}; // class ReflectionOperator

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_OPERATOR_H
