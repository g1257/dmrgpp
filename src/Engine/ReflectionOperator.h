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

template<typename RealType>
class ReflectionItem {

public:

	enum { DIAGONAL,PLUS,MINUS};

	ReflectionItem(size_t i)
	: type_(DIAGONAL),i_(i),vec_(0),value_(0)
	{}

	ReflectionItem(size_t type1,size_t i,const std::vector<RealType>& v)
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
				RealType val =vec_[k];
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
			RealType val = vec_[k];
			if (isAlmostZero(val)) continue;
			transform.pushCol(k);
			transform.pushValue(val);
			counter++;
		}
	}

	bool operator==(const ReflectionItem<RealType>& item2) const
	{
		if (type_>item2.type_) return (item2==*this);

		if (type_==DIAGONAL) {
			if (item2.type_==DIAGONAL) {
				return (i_==item2.i_);
			}
			return false;
		}

		std::vector<RealType> v3=(-1.0)*vec_;
		return (equalV(vec_,item2.vec_) || equalV(v3,item2.vec_));
	}

private:

	bool equalV(const std::vector<RealType>& v1,const std::vector<RealType>& v2) const
	{
		for (size_t i=0;i<v1.size();i++) {
			RealType x = v1[i]-v2[i];
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
	std::vector<RealType> vec_;
	RealType value_;
}; // class ReflectionItem


template<typename LeftRightSuperType>
class ReflectionOperator {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::SparseMatrixType
			SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef ReflectionItem<RealType> ItemType;
	typedef std::vector<ComplexOrRealType> VectorType;

public:

	ReflectionOperator(const LeftRightSuperType& lrs,size_t n0,bool isEnabled)
	: lrs_(lrs),
	  n0_(n0),
	  progress_("ReflectionOperator",0),
	  plusSector_(0),
	  isEnabled_(isEnabled),
	  reflected_(n0_,n0_)
	{
		size_t counter=0;
		for (size_t i=0;i<reflected_.rank();i++) {
			reflected_.setRow(i,counter);
			reflected_.pushCol(i);
			reflected_.pushValue(1);
			counter++;
		}
		reflected_.setRow(reflected_.rank(),counter);
	}

	void check(const std::vector<size_t>& sectors)
	{
		if (!isEnabled_) return;
		PsimagLite::CrsMatrix<RealType> sSuper;
		setS(sSuper);
		updateReflected(sSuper);
		PsimagLite::CrsMatrix<RealType> sSector;
		extractCurrentSector(sSector,sSuper,sectors);
		std::vector<ItemType> items;
		computeItems(items,sSector);
		std::vector<ItemType> items2;
		makeUnique(items2,items);
		assert(items2.size()==sSector.rank());
		setTransform(items2);
	}

	void changeBasis(const PsimagLite::Matrix<ComplexOrRealType>& transform)
	{
		assert(false);
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
		printFullMatrix(rT,"ConjTranform");
		SparseMatrixType matrix2;
		printFullMatrix(matrix,"OriginalHam");
		multiply(matrix2,transform_,tmp);
		printFullMatrix(transform_,"transform");
		printFullMatrix(matrix2,"FinalHam");
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

	void setS(PsimagLite::CrsMatrix<RealType>& snew)
	{
		size_t total = lrs_.super().size();
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		size_t counter = 0;
		RealType sign = 1.0;
		snew.resize(total);
		for (size_t i=0;i<total;i++) {
			snew.setRow(i,counter);
			size_t x=0,y=0;
			pack1.unpack(x,y,lrs_.super().permutation(i));

			size_t x0=0,x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			size_t y0=0,y1=0;
			pack3.unpack(y0,y1,lrs_.right().permutation(y));

			for (int k1=reflected_.getRowPtr(x0);k1<reflected_.getRowPtr(x0+1);k1++) {
				size_t x0r = reflected_.getCol(k1);
				RealType val1 = reflected_.getValue(k1);
				for (int k2=reflected_.getRowPtr(y1);k2<reflected_.getRowPtr(y1+1);k2++) {
					size_t y1r = reflected_.getCol(k2);
					RealType val2 = reflected_.getValue(k2);
					size_t xprime = pack2.pack(y1r,y0,lrs_.left().permutationInverse());
					size_t yprime = pack3.pack(x1,x0r,lrs_.right().permutationInverse());
					size_t iprime = pack1.pack(xprime,yprime,lrs_.super().permutationInverse());
			//			size_t signCounter=0;
			//			if (x0==3) signCounter++;
			//			if (x1==3) signCounter++;
			//			if (y0==3) signCounter++;
			//			if (y1==3) signCounter++;
			//			sign = 1;
			//			if (signCounter&1) sign=-1;

					snew.pushCol(iprime);
					snew.pushValue(sign*val1*val2);
					counter++;
				}
			}
		}
		snew.setRow(total,counter);
	}

	void extractCurrentSector(PsimagLite::CrsMatrix<RealType>& sSector,
				  const PsimagLite::CrsMatrix<RealType>& sSuper,
				  const std::vector<size_t>& sectors) const
	{
		assert(sectors.size()==1);
		size_t m = sectors[0];
		size_t offset = lrs_.super().partition(m);
		size_t total = lrs_.super().partition(m+1)-offset;
		sSector.resize(total);
		size_t counter = 0;
		for (size_t i=0;i<total;i++) {
			sSector.setRow(i,counter);
			for (int k=sSuper.getRowPtr(i+offset);k<sSuper.getRowPtr(i+offset+1);k++) {
				size_t col = sSuper.getCol(k);
				RealType val = sSuper.getValue(k);
				assert(col>=offset && col<total+offset);
				sSector.pushCol(col-offset);
				sSector.pushValue(val);
				counter++;
			}
		}
		sSector.setRow(total,counter);
	}

	void updateReflected(const PsimagLite::CrsMatrix<RealType>& sSuper)
	{
		// for each x0 find the reflected of
		// x0 0 0 x0
		// decompose into
		// x0' 0 0 x0'
		// then x0'=reflected(x0)
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);
//		PackIndicesType pack2(ns/n0_);
//		PackIndicesType pack3(n0_);
		reflected_.resize(ns);
		size_t counter = 0;
		for (size_t x0=0;x0<ns;x0++) {
			reflected_.setRow(x0,counter);
			//size_t x = pack2.pack(x0,0,lrs_.left().permutationInverse());
			//size_t y = pack3.pack(0,x0,lrs_.right().permutationInverse());
			size_t i = pack1.pack(x0,x0,lrs_.super().permutationInverse());
			for (int k=sSuper.getRowPtr(i);k<sSuper.getRowPtr(i+1);k++) {
				size_t col = sSuper.getCol(k);
				RealType val = sSuper.getValue(k);
				assert(!isAlmostZero(val));
				size_t xprime=0,yprime=0;
				pack1.unpack(xprime,yprime,lrs_.super().permutation(col));

				//size_t x0prime=0,x1prime=0;
				//pack2.unpack(x0prime,x1prime,lrs_.left().permutation(xprime));

				reflected_.pushCol(xprime);
				reflected_.pushValue(val);
				counter++;
			}
		}
		reflected_.setRow(ns,counter);
	}

//	size_t findReflected(size_t x) const
//	{
//		if (s_.rank()==0) return x;
//		for (size_t i=0;i<s_.rank();i++) {
//			for (int k=s_.getRowPtr(i);k<s_.getRowPtr(i+1);k++) {
//				size_t col = s_.getCol(k);
//				if (col==x) return i;
//			}
//		}
//		std::string s(__FILE__);
//		s += " findReflected\n";
//		throw std::runtime_error(s.c_str());
//	}

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

	void computeItems(std::vector<ItemType>& items,
			  PsimagLite::CrsMatrix<RealType>& sSector)
	{
		printFullMatrix(sSector,"sSector");

		for (size_t i=0;i<sSector.rank();i++) {
			std::vector<RealType> v(sSector.rank(),0);
			int equalCounter=0;
			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
				size_t col = sSector.getCol(k);
				if (isAlmostZero(sSector.getValue(k))) continue;
				v[col] =  sSector.getValue(k);
				if (col==i) equalCounter++;
				else equalCounter--;
			}
			if (equalCounter==1) {
				ItemType item(i); //v[i]);
				items.push_back(item);
				continue;
			}

			// add plus
			ItemType item(ItemType::PLUS,i,v);
			items.push_back(item);
			// add minus
			ItemType item2(ItemType::MINUS,i,v);
			items.push_back(item2);
		}
	}

	void setTransform(const std::vector<ItemType>& buffer)
	{
		transform_.resize(buffer.size());
		assert(buffer.size()==transform_.rank());
		size_t counter = 0;
		size_t row = 0;

		for (size_t i=0;i<buffer.size();i++) {
			if (buffer[i].type()==ItemType::MINUS) continue;
			transform_.setRow(row++,counter);
			buffer[i].setTransformPlus(transform_,counter);
		}

		for (size_t i=0;i<buffer.size();i++) {
			if (buffer[i].type()!=ItemType::MINUS) continue;
			transform_.setRow(row++,counter);
			buffer[i].setTransformMinus(transform_,counter);
		}
		transform_.setRow(transform_.rank(),counter);
	}

	void makeUnique(std::vector<ItemType>& dest,std::vector<ItemType>& src)
	{
		size_t zeros=0;
		size_t pluses=0;
		size_t minuses=0;
		for (size_t i=0;i<src.size();i++) {
			src[i].postFix();
		}

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
			if (item.type()==ItemType::DIAGONAL) zeros++;
			if (item.type()==ItemType::PLUS) pluses++;
			if (item.type()==ItemType::MINUS) minuses++;

			dest.push_back(item);
		}
		std::ostringstream msg;
		msg<<pluses<<" +, "<<minuses<<" -, "<<zeros<<" zeros.";
		progress_.printline(msg,std::cout);
		plusSector_ = zeros + pluses;
	}

	template<typename SomeFieldType>
	void printFullMatrix(const PsimagLite::CrsMatrix<SomeFieldType>& s,const std::string& name) const
	{
		PsimagLite::Matrix<ComplexOrRealType> fullm(s.rank(),s.rank());
		crsMatrixToFullMatrix(fullm,s);
		std::cout<<"--------->   "<<name<<" <----------\n";
		symbolicPrint(std::cout,fullm);
		//std::cout<<fullm;

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
	PsimagLite::CrsMatrix<RealType> reflected_;
	SparseMatrixType transform_;
}; // class ReflectionOperator

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_OPERATOR_H
