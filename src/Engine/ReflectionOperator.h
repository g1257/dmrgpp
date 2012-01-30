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
	  reflectionMap_(n0_)
	{
		for (size_t i=0;i<reflectionMap_.size();i++) reflectionMap_[i] = i;
	}

	void check(const std::vector<size_t>& sectors)
	{
		if (!isEnabled_) return;
		std::vector<size_t> mapL,mapR;
		setMap1(mapL,mapR);
		//reflectionMap_ = reflectionMap;
		SparseMatrixType sSector;
		setSsector(sSector,mapL,mapR,sectors);
		computeItems(sSector);
		checkTransform(sSector);
	}

	void changeBasis(const PsimagLite::Matrix<ComplexOrRealType>& transform,size_t direction)
	{
		if (!isEnabled_) return;
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

	void setMap1(std::vector<size_t>& mapL,std::vector<size_t>& mapR)
	{
		size_t ns = lrs_.left().size();
		PackIndicesType pack2(ns/n0_);
		PackIndicesType pack3(n0_);
		mapL.resize(ns);
//		std::vector<size_t> reflectionMapInverse(reflectionMap_.size());
//		for (size_t i=0;i<reflectionMapInverse.size();i++)
//			reflectionMapInverse[reflectionMap_[i]]=i;

		for (size_t i=0;i<mapL.size();i++) {
			size_t x0 = 0, x1 = 0;
			pack2.unpack(x0,x1,lrs_.left().permutation(i));
			size_t x0r = x0;
			mapL[i] = pack3.pack(x1,x0r,lrs_.right().permutationInverse());
		}

		mapR.resize(ns);
		for (size_t i=0;i<mapR.size();i++) {
			size_t x0 = 0, x1 = 0;
			pack3.unpack(x0,x1,lrs_.right().permutation(i));
			size_t x0r = x0;
			mapR[i] = pack2.pack(x1,x0r,lrs_.left().permutationInverse());
		}
	}

	void setSsector(SparseMatrixType& sSector,
			const std::vector<size_t>& mapL,
			const std::vector<size_t>& mapR,
			const std::vector<size_t>& sectors)
	{
		assert(sectors.size()==1);
		size_t m = sectors[0];
		size_t offset = lrs_.super().partition(m);
		size_t total = lrs_.super().partition(m+1)-offset;
		sSector.resize(total);
		size_t counter = 0;
		size_t ns = lrs_.left().size();
		PackIndicesType pack1(ns);

//		std::vector<size_t> reflectionMapInverse(reflectionMap.size());
//		for (size_t i=0;i<reflectionMapInverse.size();i++)
//			reflectionMapInverse[reflectionMap[i]]=i;

		for (size_t i=0;i<total;i++) {
			std::cout<<"QN for i="<<i<<" is "<<lrs_.super().qn(i+offset)<<"\n";
			sSector.setRow(i,counter);
			// i is reflected into one and only one j
			size_t x = 0, y = 0;
			pack1.unpack(x,y,lrs_.super().permutation(i+offset));
			size_t xprime = mapL[x];
			size_t yprime = mapR[y];
			size_t j = pack1.pack(yprime,xprime,lrs_.super().permutationInverse());
			assert(j>=offset && j<total+offset);
			sSector.pushValue(1);
			sSector.pushCol(j-offset);
			counter++;
		}
		sSector.setRow(total,counter);
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
	std::vector<size_t> reflectionMap_;
	SparseMatrixType transform_;
}; // class ReflectionOperator

} // namespace Dmrg 

/*@}*/
#endif // REFLECTION_OPERATOR_H
