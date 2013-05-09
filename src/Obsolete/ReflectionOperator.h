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
 *  Critical problems:
 *
 *  - support for fermionic reflections
 *
 *  Low priority work that needs to be done:
 *
 *  - support for bases that change from site to site
 *
 */
#ifndef REFLECTION_OPERATOR_H
#define REFLECTION_OPERATOR_H

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "LAPACK.h"
#include "Sort.h"
#include "ReflectionTransform.h"

namespace Dmrg {

template<typename LeftRightSuperType,typename ConcurrencyType>
class ReflectionOperator {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::SparseMatrixType
			SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef SparseVector<typename VectorType::value_type> SparseVectorType;
	typedef ReflectionTransform<RealType,SparseMatrixType> ReflectionTransformType;

	enum {AVAILABLE,NOT_AVAILABLE,COLOR};

public:

	ReflectionOperator(LeftRightSuperType& lrs,
			   ConcurrencyType& concurrency,
			   size_t n0,
			   bool isEnabled,
			   size_t expandSys)
	: lrs_(lrs),
	  concurrency_(concurrency),
	  n0_(n0),
	  progress_("ReflectionOperator",0),
	  isEnabled_(isEnabled),
	  expandSys_(expandSys),
	  reflectedLeft_(n0_,n0_),
	  reflectedRight_(n0_,n0_),
	  reflectionTransform_(false)
	{
		size_t counter=0;
		for (size_t i=0;i<reflectedLeft_.rank();i++) {
			reflectedLeft_.setRow(i,counter);
			reflectedLeft_.pushCol(i);
			reflectedLeft_.pushValue(1);
			counter++;
		}

		reflectedLeft_.setRow(reflectedLeft_.rank(),counter);
		reflectedRight_=reflectedLeft_;
	}

	void update(const typename PsimagLite::Vector<size_t>::Type& sectors)
	{
		if (!isEnabled_) return;
//		SparseMatrixType sSuper;
//		setS(sSuper);
		SparseMatrixType sSector;
		setSsector(sSector,sectors);
		updateReflected();
//		extractCurrentSector(sSector,sSuper,sectors);
		reflectionTransform_.update(sSector);
	}

	template<typename SomeStructType>
	void updateKeptStates(size_t& keptStates,
			      const SomeStructType& cacheLeft,
			      const SomeStructType& cacheRight)
	{
		const PsimagLite::Matrix<ComplexOrRealType>& transform1 = cacheLeft.transform;
		const PsimagLite::Matrix<ComplexOrRealType>& transform2 = cacheRight.transform;

		if (!isEnabled_) return;
		if (keptStates>=transform1.n_col()) return;
		PsimagLite::OstringStream msg;
		msg<<"updateKeptStates";
		progress_.printline(msg,std::cout);

		check(cacheLeft.removedIndices,reflectedLeft_,transform1,transform2);
		check(cacheRight.removedIndices,reflectedRight_,transform2,transform1);

	}

	void transform(SparseMatrixType& matrixA,
		       SparseMatrixType& matrixB,
		       const SparseMatrixType& matrix) const
	{
		assert(isEnabled_);
		reflectionTransform_.transform(matrixA,matrixB,matrix);
	}

	template<typename SomeVectorType>
	void setInitState(const SomeVectorType& initVector,
			  SomeVectorType& initVector1,
			  SomeVectorType& initVector2) const
	{
		assert(isEnabled_);
		return reflectionTransform_.setInitState(initVector,initVector1,initVector2);
	}

	RealType setGroundState(VectorType& gs,
				const RealType& gsEnergy1,
				const VectorType& gsVector1,
				const RealType& gsEnergy2,
				const VectorType& gsVector2) const
	{
		assert(isEnabled_);
		if (gsEnergy1<=gsEnergy2) {
			reflectionTransform_.setGs(gs,gsVector1,1.0);
			return gsEnergy1;
		}
		reflectionTransform_.setGs(gs,gsVector2,-1.0);
		return gsEnergy2;
	}

	const LeftRightSuperType& leftRightSuper() const { return lrs_; }

	bool isEnabled() const { return isEnabled_; }

	void changeBasis(const PsimagLite::Matrix<ComplexOrRealType>& transform1,
			 const PsimagLite::Matrix<ComplexOrRealType>& transform2)
	{
		if (!isEnabled_) return;

		SparseMatrixType newreflected;

		changeBasis(newreflected,reflectedLeft_,transform1,transform2);
//		if (newreflected.rank()!=reflectedLeft_.rank()) {
//			printFullMatrix(newreflected,"newreflectedLeft",0,1e-6);
//			transform1.print(std::cerr,1e-6);
//			transform2.print(std::cerr,1e-6);
//		}
		reflectedLeft_ = newreflected;
		//normalize(reflectedLeft_);

		changeBasis(newreflected,reflectedRight_,transform2,transform1);
		reflectedRight_ = newreflected;
		//normalize(reflectedRight_);

//		diagBasis();
	}

	void diagBasis()
	{
//		SparseMatrixType sSector;
//		setSsector(sSector);
//		reflectionTransform_.update(sSector);

//		const SparseMatrixType& Q1 = reflectionTransform_.getTransform(0);
//		const SparseMatrixType& Qm = reflectionTransform_.getTransform(1);

//		SparseMatrixType transf;
//		transformPartialLeft(transf,Q1,Qm);
//		PsimagLite::Matrix<ComplexOrRealType> fullm;
//		crsMatrixToFullMatrix(fullm,transf);
//		lrs_.leftNonConst().changeBasisDirect(fullm,concurrency_);

//		transformPartialRight(transf,Q1,Qm);
//		crsMatrixToFullMatrix(fullm,transf);
//		lrs_.rightNonConst().changeBasisDirect(fullm,concurrency_);

	}

private:

//	void transformPartialLeft(SparseMatrixType& tr,const SparseMatrixType& Q1,const SparseMatrix& Qm) const
//	{

//		size_t ns = lrs_.left().size();
//		PackIndicesType pack2(ns/n0_);
//		PackIndicesType pack3(n0_);
//		PackIndicesType pack1(ns);
//		assert(reflectedLeft_.rank()==ns/n0_);
//		assert(reflectedRight_.rank()==ns/n0_);
//		typename PsimagLite::Vector<int>::Type ptr(total,-1);
//		typename PsimagLite::Vector<size_t>::Type index(total,0);
//		typename PsimagLite::Vector<ComplexOrRealType>::Type temp(total,0);

//		size_t counter = 0;
//		for (size_t i=0;i<total;i++) {
//			size_t x = 0, y = 0;
//			pack1.unpack(x,y,lrs_.super().permutation(i));
//			tr.setRow(x,counter);
//			for (int k=Q1.getRowPtr(i);k<Q1.getRowPtr(i+1);k++) {
//				size_t col = Q1.getCol(k);
//				pack1.unpack(xprime,yprime,lrs_.super().permutation(i));


//	}

//	void transformPartial(SparseMatrixType& s1,const PsimagLite::Matrix<ComplexOrRealType>& fullm)
//	{
//		SparseMatrixType m1(fullm);
//		SparseMatrixType m1Conj;
//		transposeConjugate(m1Conj,m1);
//		SparseMatrixType tmp = s1*m1Conj;
//		s1 = m1*tmp;
//		s1.checkValidity();
//	}

//	void transformPartial(SparseMatrixType& s1,const typename PsimagLite::Vector<RealType>::Type& eigs)
//	{
//		size_t n = s1.rank();
//		s1.resize(n);
//		size_t counter = 0;
//		for (size_t i=0;i<n;i++) {
//			s1.setRow(i,counter);
//			s1.pushCol(i);
//			s1.pushValue(eigs[i]);
//			counter++;
//		}
//		s1.setRow(n,counter);
//		s1.checkValidity();
//	}

	void check(const typename PsimagLite::Vector<size_t>::Type& removedIndices,
		   const SparseMatrixType& reflected,
		   const PsimagLite::Matrix<ComplexOrRealType>& transform1,
		   const PsimagLite::Matrix<ComplexOrRealType>& transform2)
	{

		SparseMatrixType newreflected;

		changeBasis(newreflected,reflected,transform1,transform2);

//		RealType eps = 1e-6;
//		printFullMatrix(newreflected,"newreflected",0,eps);


		typename PsimagLite::Vector<size_t>::Type x;
		for (size_t ii=0;ii<removedIndices.size();ii++) {
			size_t i = removedIndices[ii];
			for (int k=newreflected.getRowPtr(i);k<newreflected.getRowPtr(i+1);k++) {
				ComplexOrRealType val = newreflected.getValue(k);
				if (isAlmostZero(val,1e-8)) continue;
				size_t col = newreflected.getCol(k);
				x.push_back(col);
			}
		}
		Sort<typename PsimagLite::Vector<size_t>::Type > sort;
		typename PsimagLite::Vector<size_t>::Type iperm(x.size());
		sort.sort(x,iperm);

		typename PsimagLite::Vector<size_t>::Type diffs;
		getDifferences(diffs,x,removedIndices);
		for (size_t i=0;i<diffs.size();i++)
			std::cerr<<"diffs["<<i<<"]="<<diffs[i]<<"\n";

		for (size_t i=0;i<x.size();i++)
			std::cerr<<"x["<<i<<"]="<<x[i]<<"\n";

		for (size_t i=0;i<removedIndices.size();i++)
			std::cerr<<"removed["<<i<<"]="<<removedIndices[i]<<"\n";

	}

	void getDifferences(typename PsimagLite::Vector<size_t>::Type& diffs,
			    const typename PsimagLite::Vector<size_t>::Type& x1,
			    const typename PsimagLite::Vector<size_t>::Type& x2) const
	{
		typename PsimagLite::Vector<size_t>::const_iterator::Type it2 = x2.begin();
		for (size_t i=0;i<x1.size();i++) {
			typename PsimagLite::Vector<size_t>::const_iterator::Type it = find(it2,x2.end(),x1[i]);
			if (it == x2.end()) {
				diffs.push_back(x1[i]);
				continue;
			}
			it2 = it+1;
		}
	}

	void normalize(SparseMatrixType& A) const
	{
		size_t n = A.rank();
		typename PsimagLite::Vector<ComplexOrRealType>::Type sum(n,0.0);
		for (size_t i=0;i<n;i++) {
			for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				ComplexOrRealType val = A.getValue(k);
				sum[i] += std::conj(val)*val;
			}
		}

		for (size_t i=0;i<n;i++) {
			//assert(isAlmostZero(sum[i]-1.0,1e-6));
			assert(!isAlmostZero(sum[i],1e-6));
			ComplexOrRealType x = 1.0/sqrt(sum[i]);
			for (int k=A.getRowPtr(i);k<A.getRowPtr(i+1);k++)
				A.setValues(k,A.getValue(k)*x);
		}
	}

	void changeBasis(SparseMatrixType& newreflected,
			 const SparseMatrixType& reflected,
			 const PsimagLite::Matrix<ComplexOrRealType>& transform1,
			 const PsimagLite::Matrix<ComplexOrRealType>& transform2)
	{
		assert(reflected.rank()==transform1.n_row());
		assert(reflected.rank()==transform2.n_row());

		size_t total = transform1.n_col();
		newreflected.resize(total);
		typename PsimagLite::Vector<int>::Type ptr(total,-1);
		typename PsimagLite::Vector<size_t>::Type index(total,0);
		typename PsimagLite::Vector<ComplexOrRealType>::Type temp(total,0);
		std::cerr<<"transform1="<<transform1.n_row()<<"x"<<transform1.n_col()<<"\n";
		std::cerr<<"transform2="<<transform2.n_row()<<"x"<<transform2.n_col()<<"\n";

		size_t counter = 0;
		for (size_t x=0;x<total;x++) {
			newreflected.setRow(x,counter);

			size_t itemp = 0;
			for (size_t xprime=0;xprime<transform1.n_row();xprime++) {
				ComplexOrRealType wl1 =  transform1(xprime,x);
//				ComplexOrRealType wl1 =  transform(x,xprime);
				for (int k=reflected.getRowPtr(xprime);k<reflected.getRowPtr(xprime+1);k++) {
					size_t xsecond = reflected.getCol(k);
					ComplexOrRealType r = reflected.getValue(k);
					for (size_t xthird=0;xthird<transform2.n_col();xthird++) {
//						ComplexOrRealType val = wl1 * r * transform(xthird,xsecond);
						ComplexOrRealType val = wl1 * r * transform2(xsecond,xthird);
						//if (isAlmostZero(val)) continue;
						if (ptr[xthird]<0) {
							ptr[xthird] = itemp;
							temp[ptr[xthird]] = val;
							index[ptr[xthird]] = xthird;
							itemp++;
						} else {
							temp[ptr[xthird]] += val;
						}
					}
				}
			}
			for (size_t s=0;s<itemp;s++) {
				newreflected.pushCol(index[s]);
				newreflected.pushValue(temp[s]);
				ptr[index[s]] = -1;
			}
			counter += itemp;
		}
		newreflected.setRow(newreflected.rank(),counter);
		newreflected.checkValidity();
	}

	void setSsector(SparseMatrixType& sSector,const typename PsimagLite::Vector<size_t>::Type& sectors) const
	{
		assert(sectors.size()==1);
		size_t m = sectors[0];
		size_t offset = lrs_.super().partition(m);
		size_t total = lrs_.super().partition(m+1)-offset;
		setSsector(sSector,total,offset);
	}

	void setSsector(SparseMatrixType& sSector,size_t total=0,size_t offset=0) const
	{
		if (total==0) total=lrs_.super().size();
		sSector.resize(total);
		size_t counter = 0;
		size_t ns = lrs_.left().size();
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		PackIndicesType pack1(ns);
		assert(reflectedLeft_.rank()==ns/n0_);
		assert(reflectedRight_.rank()==ns/n0_);
		typename PsimagLite::Vector<int>::Type ptr(total,-1);
		typename PsimagLite::Vector<size_t>::Type index(total,0);
		typename PsimagLite::Vector<ComplexOrRealType>::Type temp(total,0);

		for (size_t i=0;i<total;i++) {
			sSector.setRow(i,counter);
			size_t x = 0, y = 0;
			pack1.unpack(x,y,lrs_.super().permutation(i+offset));

			size_t x0=0,x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			size_t y0=0,y1=0;
			pack3.unpack(y0,y1,lrs_.right().permutation(y));

			size_t itemp = 0;

			for (int k=reflectedLeft_.getRowPtr(x0);k<reflectedLeft_.getRowPtr(x0+1);k++) {
				ComplexOrRealType val1 = reflectedLeft_.getValue(k);
				for (int k2=reflectedRight_.getRowPtr(y1);k2<reflectedRight_.getRowPtr(y1+1);k2++) {
					ComplexOrRealType val2 = val1 * reflectedRight_.getValue(k2);
					if (isAlmostZero(val2)) continue;
					size_t x0prime = reflectedLeft_.getCol(k);
					size_t xprime = pack3.pack(x1,x0prime,lrs_.right().permutationInverse());

					size_t y1prime = reflectedRight_.getCol(k2);
					size_t yprime = pack2.pack(y1prime,y0,lrs_.left().permutationInverse());

					size_t iprime = pack1.pack(yprime,xprime,lrs_.super().permutationInverse());
					assert(iprime>=offset && iprime<offset+total);
					size_t col = iprime - offset;
					if (ptr[col]<0) {
						ptr[col] = itemp;
						temp[ptr[col]] = val2;
						index[ptr[col]] = col;
						itemp++;
					} else {
						temp[ptr[col]] += val2;
					}
				}
			}
			ComplexOrRealType val = 0.0;
			for (size_t s=0;s<itemp;s++) {
				sSector.pushCol(index[s]);
				sSector.pushValue(temp[s]);
				val += std::conj(temp[s])*temp[s];
				ptr[index[s]] = -1;
				counter++;
			}
			if (isAlmostZero(val,1e-8)) {
				std::cerr<<"i="<<i<<" x0="<<x0<<" x1="<<x1<<" y0="<<y0<<" y1="<<y1<<"\n";
			}
		}
		sSector.setRow(total,counter);
		sSector.checkValidity();

//		printSparseMatrix(sSector,"sSector",1e-3);
#ifndef NDEBUG
		SparseMatrixType tmp;
		multiply(tmp,sSector,sSector);
		bool b = isTheIdentity(tmp,1e-5);
		if (b) return;
		hasZeroRows(sSector,true);
		printFullMatrix(sSector,"sSector");
		printFullMatrix(tmp,"sSector*sSector");
		assert(false);
#endif
	}

	bool hasZeroRows(const SparseMatrixType& sSector,bool verbose) const
	{
		bool flag = false;
		size_t n = sSector.rank();
		for (size_t i=0;i<n;i++) {
			ComplexOrRealType sum = 0.0;
			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
				ComplexOrRealType val = sSector.getValue(k);
				sum += std::conj(val)*val;
			}
			if (isAlmostZero(sum,1e-4)) {
				flag=true;
				if (verbose) std::cerr<<"zero row="<<i<<"\n";
			}
		}
		return flag;
	}


	void printSparseMatrix(SparseMatrixType& s1,const PsimagLite::String& label,const RealType& eps) const
	{
		std::cout<<label<<"\n";
		for (size_t i=0;i<s1.rank();i++) {
			for (int k=s1.getRowPtr(i);k<s1.getRowPtr(i+1);k++) {
				ComplexOrRealType val = s1.getValue(k);
				if (isAlmostZero(val,eps)) continue;
				std::cout<<(i+1)<<" "<<(1+s1.getCol(k))<<" "<<val<<"\n";
			}
		}
	}

	void truncate(SparseMatrixType& reflectedFinal,
		      const SparseMatrixType& reflected,
		      size_t keptstates)
	{
//		size_t n = reflected.rank();
//		if (keptstates >= n) {
			reflectedFinal = reflected;
			return;
//		}

//		reflectedFinal.resize(keptstates);
//		size_t counterl = 0;

//		typename PsimagLite::Vector<ComplexOrRealType>::Type sum(keptstates,0.0);
//		for (size_t i=0;i<keptstates;i++) {
//			reflectedFinal.setRow(i,counterl);
//			for (int k = reflected.getRowPtr(i); k < reflected.getRowPtr(i+1);k++) {
//				size_t col = reflected.getCol(k);
//				if (col>=keptstates) continue;
//				ComplexOrRealType val = reflected.getValue(k);
//				reflectedFinal.pushCol(col);
//				reflectedFinal.pushValue(val);
//				counterl++;
//				sum[i] += std::conj(val) * val;
//			}

//		}
//		reflectedFinal.setRow(reflectedFinal.rank(),counterl);
//		reflectedFinal.checkValidity();

//		// normalize
//		for (size_t i=0;i<reflectedFinal.rank();i++) {
//			assert(!isAlmostZero(sum[i],1e-8));
//			ComplexOrRealType x = 1.0/sqrt(sum[i]);
//			for (int k = reflectedFinal.getRowPtr(i); k < reflectedFinal.getRowPtr(i+1);k++)
//				reflectedFinal.setValues(k,reflectedFinal.getValue(k)*x);
//		}
//		reflectedFinal.checkValidity();
	}

	void updateReflected()
	{
		size_t ns = lrs_.left().size();
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		SparseMatrixType reflectedLeft(ns,ns);

		size_t counter = 0;
		for (size_t x=0;x<ns;x++) {
			reflectedLeft.setRow(x,counter);
			size_t x0 = 0, x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));
			for (int k=reflectedLeft_.getRowPtr(x0);k<reflectedLeft_.getRowPtr(x0+1);k++) {
				size_t x0r = reflectedLeft_.getCol(k);
				size_t col = pack3.pack(x1,x0r,lrs_.right().permutationInverse());
				ComplexOrRealType val = reflectedLeft_.getValue(k);
				//if (isAlmostZero(val)) continue;
				reflectedLeft.pushCol(col);
				reflectedLeft.pushValue(val);
				counter++;
			}
		}
		reflectedLeft.setRow(ns,counter);
		reflectedLeft.checkValidity();
		reflectedLeft_ = reflectedLeft;

		SparseMatrixType reflectedRight(ns,ns);
		counter=0;
		for (size_t x=0;x<ns;x++) {
			size_t x0 = 0, x1=0;
			reflectedRight.setRow(x,counter);
			pack3.unpack(x0,x1,lrs_.right().permutation(x));
			for (int k=reflectedRight_.getRowPtr(x1);k<reflectedRight_.getRowPtr(x1+1);k++) {
				size_t x1r = reflectedRight_.getCol(k);
				ComplexOrRealType val = reflectedRight_.getValue(k);
				//if (isAlmostZero(val)) continue;
				size_t col = pack2.pack(x1r,x0,lrs_.left().permutationInverse());
				reflectedRight.pushCol(col);
				reflectedRight.pushValue(val);
				counter++;
			}
		}
		reflectedRight.setRow(ns,counter);
		reflectedRight.checkValidity();
		reflectedRight_ = reflectedRight;
	}

	LeftRightSuperType& lrs_;
	ConcurrencyType& concurrency_;
	size_t n0_; // hilbert size of one site
	PsimagLite::ProgressIndicator progress_;
	bool isEnabled_;
	size_t expandSys_;
	SparseMatrixType reflectedLeft_,reflectedRight_;
	ReflectionTransformType reflectionTransform_;
}; // class ReflectionOperator

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_OPERATOR_H
