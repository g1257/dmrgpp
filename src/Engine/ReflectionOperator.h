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
 *  - getting a l.i. set in an efficient way
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
#include "LinearlyIndependentSet.h"
#include "LAPACK.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class ReflectionOperator {

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename LeftRightSuperType::SparseMatrixType
			SparseMatrixType;
	typedef typename LeftRightSuperType::RealType RealType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
//	typedef ReflectionItem<RealType,ComplexOrRealType> ItemType;
	typedef std::vector<ComplexOrRealType> VectorType;

public:

	ReflectionOperator(const LeftRightSuperType& lrs,size_t n0,bool isEnabled,size_t expandSys)
	: lrs_(lrs),
	  n0_(n0),
	  progress_("ReflectionOperator",0),
	  plusSector_(0),
	  isEnabled_(isEnabled),
	  expandSys_(expandSys),
	  reflectedLeft_(n0_,n0_),
	  reflectedRight_(n0_,n0_),
	  seMap_(n0_)
	{
		size_t counter=0;
		for (size_t i=0;i<reflectedLeft_.rank();i++) {
			reflectedLeft_.setRow(i,counter);
			reflectedLeft_.pushCol(i);
			reflectedLeft_.pushValue(1);
			counter++;
		}
		for (size_t i=0;i<seMap_.size();i++) seMap_[i] = i;

		reflectedLeft_.setRow(reflectedLeft_.rank(),counter);
		reflectedRight_=reflectedLeft_;
	}

	void check(const std::vector<size_t>& sectors)
	{
		if (!isEnabled_) return;
		SparseMatrixType sSuper;
		setS(sSuper);
		updateReflected(sSuper);
		SparseMatrixType sSector;
		extractCurrentSector(sSector,sSuper,sectors);
		computeItems(sSector);
		checkTransform(sSector);
	}

	void changeBasis(const PsimagLite::Matrix<ComplexOrRealType>& transform,size_t direction)
	{
		if (!isEnabled_) return;
		SparseMatrixType newreflected;
		if (direction==expandSys_) {
			changeBasis(newreflected,reflectedLeft_,transform);
			reflectedLeft_ = newreflected;
		} else {
			changeBasis(newreflected,reflectedRight_,transform);
			reflectedRight_ = newreflected;
		}
		SparseMatrixType one1;
		multiply(one1,reflectedLeft_,reflectedLeft_);
		printFullMatrix(one1,"Should be I for ReflectedLeft:");
		assert(isTheIdentity(one1));

		multiply(one1,reflectedRight_,reflectedRight_);
		printFullMatrix(one1,"Should be I for ReflectedRight:");
		assert(isTheIdentity(one1));
	}

	void transform(SparseMatrixType& matrixA,
		       SparseMatrixType& matrixB,
		       const SparseMatrixType& matrix) const
	 {
		assert(isEnabled_);

		SparseMatrixType rT;
//		invert(rT,transform_);
		transposeConjugate(rT,transform_);

		SparseMatrixType tmp;
		multiply(tmp,matrix,rT);

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

	const LeftRightSuperType& leftRightSuper() const { return lrs_; }

	bool isEnabled() const { return isEnabled_; }

private:

	void checkTransform(const SparseMatrixType& sSector)
	{
		SparseMatrixType rT;
//		invert(rT,transform_);
		transposeConjugate(rT,transform_);
		SparseMatrixType tmp3;
		multiply(tmp3,transform_,rT);
		printFullMatrix(sSector,"Ssector");

		printFullMatrix(rT,"Transform");

		SparseMatrixType tmp4;
		multiply(tmp3,sSector,rT);
		multiply(tmp4,transform_,tmp3);
		printFullMatrix(tmp4,"R S R^\\dagger");
	}

//	void invert(SparseMatrixType& dest,const SparseMatrixType& src) const
//	{
//		PsimagLite::Matrix<ComplexOrRealType> fullM;
//		crsMatrixToFullMatrix(fullM,src);
//		int m = fullM.n_row();
//		int n = m;
//		int lda = m;
//		int info = 0;
//		std::vector<int> ipiv(m);
//		psimag::LAPACK::GETRF(m,n,&fullM(0,0),lda,&(ipiv[0]),info);
//		assert(info==0);
//		n = fullM.n_row();
//		lda = m;
//		// query:
//		int lwork = -1;
//		std::vector<ComplexOrRealType> work(3);
//		psimag::LAPACK::GETRI(n,&fullM(0,0),lda,&(ipiv[0]),&(work[0]),lwork,info);
//		lwork = std::real(work[0]);
//		// actual work:
//		work.resize(lwork+5);
//		psimag::LAPACK::GETRI(n,&fullM(0,0),lda,&(ipiv[0]),&(work[0]),lwork,info);
//		assert(info==0);
//		fullMatrixToCrsMatrix(dest,fullM);
//	}

	void changeBasis(SparseMatrixType& newreflected,const SparseMatrixType& reflected,const PsimagLite::Matrix<ComplexOrRealType>& transform)
	{
		newreflected.resize(reflected.rank());
		size_t counter = 0;
		assert(reflected.rank()==transform.n_row());
		assert(reflected.rank()==transform.n_row());

		size_t total = reflected.rank();
		std::vector<int> ptr(total,-1);
		std::vector<size_t> index(total,0);
		std::vector<ComplexOrRealType> temp(total,0);

		for (size_t x=0;x<total;x++) {
			newreflected.setRow(x,counter);

			size_t itemp = 0;
			for (size_t xprime=0;xprime<transform.n_row();xprime++) {
				ComplexOrRealType wl1 =  transform(xprime,x);
//				ComplexOrRealType wl1 =  transform(x,xprime);
				for (int k=reflected.getRowPtr(xprime);k<reflected.getRowPtr(xprime+1);k++) {
					size_t xsecond = reflected.getCol(k);
					ComplexOrRealType r = reflected.getValue(k);
					for (size_t xthird=0;xthird<transform.n_col();xthird++) {
//						ComplexOrRealType val = wl1 * r * transform(xthird,xsecond);
						ComplexOrRealType val = wl1 * r * transform(xsecond,xthird);
						if (isAlmostZero(val)) continue;
						if (ptr[xthird]<0) {
							ptr[xthird] = itemp;
							temp[ptr[xthird]] = val;
							index[ptr[xthird]] = xthird;
							itemp++;
						} else {
							temp[ptr[xthird]] += val;
						}
//						newreflected.pushCol(xthird);
//						newreflected.pushValue(val);
//						counter++;
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
		newreflected.setRow(reflected.rank(),counter);
		newreflected.checkValidity();
	}

	void setSNew(SparseMatrixType& snew) const
	{
		size_t total = lrs_.super().size();
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		size_t counter = 0;
//		RealType sign = 1.0;
		snew.resize(total);
//		const SparseMatrixType& aMatrix = reflectedLeft_;
//		const SparseMatrixType& bMatrix = reflectedRight_;

		std::vector<size_t> seMap2(seMap_.size());
		for (size_t i=0;i<seMap2.size();i++)
			seMap2[seMap_[i]]=i;

		std::vector<int> ptr(total,-1);
		std::vector<size_t> index(total,0);
		std::vector<ComplexOrRealType> temp(total,0);

		for (size_t i=0;i<total;i++) {
			snew.setRow(i,counter);
			size_t x=0,y=0;
			pack1.unpack(x,y,lrs_.super().permutation(i));
			size_t col = pack1.pack(y,x,lrs_.super().permutationInverse());
			snew.pushCol(col);
			snew.pushValue(1);
			counter++;
		}
		snew.setRow(total,counter);
	}

	void setS(SparseMatrixType& snew) const
	{
		size_t total = lrs_.super().size();
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		size_t counter = 0;
		RealType sign = 1.0;
		snew.resize(total);
		const SparseMatrixType& aMatrix = reflectedLeft_;
		const SparseMatrixType& bMatrix = reflectedRight_;

		std::vector<size_t> seMap2(seMap_.size());
		for (size_t i=0;i<seMap2.size();i++)
			seMap2[seMap_[i]]=i;

		std::vector<int> ptr(total,-1);
		std::vector<size_t> index(total,0);
		std::vector<ComplexOrRealType> temp(total,0);

		for (size_t i=0;i<total;i++) {
			snew.setRow(i,counter);
			size_t x=0,y=0;
			pack1.unpack(x,y,lrs_.super().permutation(i));

			size_t x0=0,x1=0;
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			size_t y0=0,y1=0;
			pack3.unpack(y0,y1,lrs_.right().permutation(y));

			size_t itemp = 0;
			for (int k1=aMatrix.getRowPtr(x0);k1<aMatrix.getRowPtr(x0+1);k1++) {
				size_t x0r = seMap_[aMatrix.getCol(k1)];
//				size_t x0r = aMatrix.getCol(k1);
				ComplexOrRealType val1 = aMatrix.getValue(k1);
				for (int k2=bMatrix.getRowPtr(y1);k2<bMatrix.getRowPtr(y1+1);k2++) {
					size_t y1r = seMap2[bMatrix.getCol(k2)];
//					size_t y1r = bMatrix.getCol(k2);
					ComplexOrRealType val2 = sign*val1*bMatrix.getValue(k2);
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

					if (ptr[iprime]<0) {
						ptr[iprime] = itemp;
						temp[ptr[iprime]] = val2;
						index[ptr[iprime]] = iprime;
						itemp++;
					} else {
						temp[ptr[iprime]] += val2;
					}
				}
			}
			for (size_t s=0;s<itemp;s++) {
				snew.pushCol(index[s]);
				snew.pushValue(temp[s]);
				ptr[index[s]] = -1;
			}
			counter += itemp;
		}
		snew.setRow(total,counter);

		snew.checkValidity();

		SparseMatrixType one1;
		multiply(one1,snew,snew);
		printFullMatrix(one1,"Should be I for FR:");
		assert(isTheIdentity(one1));

//		PsimagLite::Matrix<ComplexOrRealType> m1;
//		crsMatrixToFullMatrix(m1,snew);
//		PsimagLite::Matrix<ComplexOrRealType> m2(m1.n_row(),m1.n_col());
//		for (size_t i=0;i<m1.n_row();i++) {
//			for (size_t j=0;j<m1.n_col();j++) {
//				m2(i,j) = 0;
//				for (size_t k=0;k<m1.n_col();k++) {
//					 m2(i,j) += m1(i,k) * m1(k,j);
//				}
//			}
//		}
//		std::cout<<"\n";
//		std::cout<<m2<<"\n";

//		SparseMatrixType tmp;
//		multiply(tmp,snew,snew);

	}

	void extractCurrentSector(SparseMatrixType& sSector,
				  const SparseMatrixType& sSuper,
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
				ComplexOrRealType val = sSuper.getValue(k);

				assert(col>=offset && col<total+offset);
				sSector.pushCol(col-offset);
				sSector.pushValue(val);
				counter++;
			}
		}
		sSector.setRow(total,counter);
		sSector.checkValidity();

		SparseMatrixType tmp;
		multiply(tmp,sSector,sSector);
		assert(isTheIdentity(tmp));
//		printFullMatrix(tmp,"ShouldBeI");

//		PsimagLite::Matrix<ComplexOrRealType> m1;
//		crsMatrixToFullMatrix(m1,sSector);
//		PsimagLite::Matrix<ComplexOrRealType> m2(m1.n_row(),m1.n_col());
//		for (size_t i=0;i<m1.n_row();i++) {
//			for (size_t j=0;j<m1.n_col();j++) {
//				m2(i,j) = 0;
//				for (size_t k=0;k<m1.n_col();k++) {
//					 m2(i,j) += m1(i,k) * m1(k,j);
//				}
//			}
//		}
//		std::cout<<"Should be I:\n";
//		std::cout<<m2<<"\n";
	}

	void updateReflected(const SparseMatrixType& sSuper)
	{
		// for each x0 find the reflected of
		// x0 0 0 x0
		// decompose into
		// x0' 0 0 x0'
		// then x0'=reflected(x0)
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);
		size_t ne = lrs_.super().size()/ns;
		assert(ns==ne);
		reflectedLeft_.resize(ns);
		reflectedRight_.resize(ne);
		size_t counter = 0;
		seMap_.resize(ns);
		for (size_t i=0;i<seMap_.size();i++) seMap_[i] = i;
		for (size_t x0=0;x0<ns;x0++) {
			reflectedLeft_.setRow(x0,counter);
			reflectedRight_.setRow(x0,counter);
			size_t i = pack1.pack(x0,x0,lrs_.super().permutationInverse());
			for (int k=sSuper.getRowPtr(i);k<sSuper.getRowPtr(i+1);k++) {
				size_t col = sSuper.getCol(k);
				ComplexOrRealType val = sSuper.getValue(k);
				if (isAlmostZero(val)) continue;
				ComplexOrRealType val1 = 0, val2 = 0;
				assert(isAlmostZero(std::real(val)-1.0));
//				if (std::real(val)<0) {
//					val1 = sqrt(-val);
//					val2 = -val1;
//				} else {
					val1 = val2 = sqrt(val);
//				}

				size_t xprime=0,yprime=0;
				pack1.unpack(xprime,yprime,lrs_.super().permutation(col));
				seMap_[xprime] = yprime;
				//assert(xprime==yprime);
				reflectedLeft_.pushCol(xprime);
				reflectedLeft_.pushValue(1);

				reflectedRight_.pushCol(yprime);
				reflectedRight_.pushValue(1);
				counter++;
			}
		}
		reflectedLeft_.setRow(ns,counter);
		reflectedRight_.setRow(ns,counter);
		reflectedLeft_.checkValidity();
		reflectedRight_.checkValidity();

//		SparseMatrixType one1;
//		multiply(one1,reflectedLeft_,reflectedLeft_);
//		printFullMatrix(one1,"Should be I for ReflectedLeft:");
//		assert(isTheIdentity(one1));

//		multiply(one1,reflectedRight_,reflectedRight_);
//		printFullMatrix(one1,"Should be I for ReflectedRight:");
//		assert(isTheIdentity(one1));


	}

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

	void computeItems(const SparseMatrixType& sSector)
	{
		printFullMatrix(sSector,"sSector");
		LinearlyIndependentSet<RealType,VectorType> lis(sSector.rank());

		for (size_t i=0;i<sSector.rank();i++) {
			std::vector<ComplexOrRealType> v(sSector.rank(),0);
			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
				size_t col = sSector.getCol(k);
				if (isAlmostZero(sSector.getValue(k))) continue;
				v[col] =  sSector.getValue(k);
			}
			v[i] += 1.0;
			lis.push(v);
		}
		plusSector_=lis.size();
		for (size_t i=0;i<sSector.rank();i++) {
			std::vector<ComplexOrRealType> v(sSector.rank(),0);
			for (int k=sSector.getRowPtr(i);k<sSector.getRowPtr(i+1);k++) {
				size_t col = sSector.getCol(k);
				if (isAlmostZero(sSector.getValue(k))) continue;
				v[col] =  sSector.getValue(k);
			}
			v[i] -= 1.0;
			lis.push(v);
		}
		size_t minuses = lis.size()-plusSector_;
		std::ostringstream msg;
		msg<<plusSector_<<" +, "<<minuses<<" -.";
		progress_.printline(msg,std::cout);
		assert(lis.size()==sSector.rank());
		lis.fill(transform_);
	}

	void printFullMatrix(const SparseMatrixType& s,const std::string& name) const
	{
		PsimagLite::Matrix<ComplexOrRealType> fullm(s.rank(),s.rank());
		crsMatrixToFullMatrix(fullm,s);
		std::cout<<"--------->   "<<name<<" <----------\n";
		try {
//			mathematicaPrint(std::cout,fullm);
//			symbolicPrint(std::cout,fullm);
		} catch (std::exception& e) {
			std::cout<<fullm;
		}

		std::cout<<fullm;

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
	size_t expandSys_;
	SparseMatrixType reflectedLeft_,reflectedRight_;
	std::vector<size_t> seMap_;
	SparseMatrixType transform_;
}; // class ReflectionOperator

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_OPERATOR_H
