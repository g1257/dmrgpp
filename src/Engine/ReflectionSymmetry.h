
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef REFLECTION_SYMM_H
#define REFLECTION_SYMM_H
#include <iostream>
#include "ProgressIndicator.h"
#include "CrsMatrix.h"
#include "Vector.h"

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

	template<typename GeometryType,typename BasisType>
	class ReflectionSymmetry  {

		typedef typename GeometryType::RealType RealType;
		typedef typename BasisType::WordType WordType;
		typedef PsimagLite::CrsMatrix<RealType> SparseMatrixType;
		typedef std::vector<RealType> VectorType;

		typedef ReflectionItem ItemType;

	public:

		ReflectionSymmetry(const BasisType& basis,const GeometryType& geometry)
		: progress_("ReflectionSymmetry",0),
		  transform_(basis.size(),basis.size()),
		  plusSector_(0)
		{
			size_t hilbert = basis.size();
			size_t numberOfDofs = basis.dofs();
			size_t numberOfSites = geometry.numberOfSites();
			size_t termId = 0;
//			size_t counter=0;
			std::vector<ItemType> buffer;
			for (size_t ispace=0;ispace<hilbert;ispace++) {
				std::vector<WordType> y(numberOfDofs,0);
				for (size_t dof=0;dof<numberOfDofs;dof++) {
					WordType x = basis(ispace,dof);
					for (size_t site=0;site<numberOfSites;site++) {
						size_t reflectedSite = geometry.findReflection(site,termId);
						size_t thisSiteContent = x & 1;
						x >>=1; // go to next site
						addTo(y[dof],thisSiteContent,reflectedSite);
						if (!x) break;
					}
				}

				size_t yIndex = basis.perfectIndex(y);
//				s_.setRow(ispace,counter);
//				s_.pushCol(yIndex);
//				s_.pushValue(1.0);
//				counter++;
				if (yIndex==ispace) { // then S|psi> = |psi>
					ItemType item1(ispace);
					buffer.push_back(item1);
					continue;
				}
				// S|psi> != |psi>
				// Add normalized +
				ItemType item2(ispace,yIndex,ItemType::PLUS);
				buffer.push_back(item2);

				// Add normalized -
				ItemType item3(ispace,yIndex,ItemType::MINUS);
				buffer.push_back(item3);
			}
//			s_.setRow(s_.rank(),counter);
			setTransform(buffer);
//			checkTransform();
		}

		void transform(SparseMatrixType& matrixA,
			       SparseMatrixType& matrixB,
			       const SparseMatrixType& matrix) const
		{
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);

			SparseMatrixType tmp;
			multiply(tmp,matrix,rT);

			SparseMatrixType matrix2;
			multiply(matrix2,transform_,tmp);

			split(matrixA,matrixB,matrix2);
		}

		RealType setGroundState(VectorType& gs,
					const RealType& gsEnergy1,
					const VectorType& gsVector1,
					const RealType& gsEnergy2,
					const VectorType& gsVector2)
		{
			size_t rank = gsVector1.size() + gsVector2.size();
			if (gsEnergy1<=gsEnergy2) {
				setGs(gs,gsVector1,rank,0);
				return gsEnergy1;
			}
			setGs(gs,gsVector2,rank,gsVector1.size());
			return gsEnergy2;
		}

	private:

		void calcReflectionOpSystem()
		{
			for (ArrangementsType i(levels,sizes);!i.end();i.next()) {
				for (ArrangementsType iprime(levels,sizes);!i.end();i.next()) {
					for (ArrangementsType mu(levels,n0);!i.end();i.next()) {
						size_t firstMu = mu[0];
						size_t lastI = i[levels-1];
						size_t lastIprime = iprime[levels-1];
						size_t lastMu = mu[levels-1];
						size_t p = permutationInverse[levels-1](lastI,lastMu);
						size_t pprime = permutationInverse[levels-1](lastIprime,firstMu);
						s_(p,pprime) += calcReflectionOpSystem(i,iprime,mu);
					}
				}
			}
		}

		void calcReflectionOpSystem(const ArrangementsType& i,
					    const ArrangementsType& iprime,
					    const ArrangementsType& mu)
		{
			for (size_t x=0;x<i.size();x++) {
				p = permutationInverse[x](i[x-1],mu[x]);
				tmp *= w[x](p,i[x]);
				pprime = permutationInverse[levels-x-1](iprime[x-1],mu[levels-x-1]);
				tmp *= w[x](pprime,iprime[x]);
			}
		}

		void setGs(VectorType& gs,const VectorType& v,size_t rank,size_t offset)
		{
			std::vector<RealType> gstmp(rank,0);

			for (size_t i=0;i<v.size();i++) {
				gstmp[i+offset]=v[i];
			}
			SparseMatrixType rT;
			transposeConjugate(rT,transform_);
			multiply(gs,rT,gstmp);
		}

		void checkTransform() const
		{
			SparseMatrixType transformTc;
			transposeConjugate(transformTc,transform_);

//			SparseMatrixType tmp;
//			multiply(tmp,s_,s_);
			PsimagLite::Matrix<RealType> mtmp;
			crsMatrixToFullMatrix(mtmp,transform_);
			std::cerr<<"&&&&&&&&&&&&&&&&&&&&&&&\n";
			std::cerr<<mtmp;
//			throw std::runtime_error("checking\n");
		}

		void addTo(WordType& yy,size_t what,size_t site) const
		{
			if (what==0) return;
			WordType mask = (1<<site);
			yy |= mask;
		}

		void setTransform(const std::vector<ItemType>& buffer2)
		{
			std::vector<ItemType> buffer;
			makeUnique(buffer,buffer2);
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

		void isIdentity(const SparseMatrixType& s,const std::string& label) const
		{
			std::cerr<<"Checking label="<<label<<"\n";
			for (size_t i=0;i<s.rank();i++) {
				for (int k=s.getRowPtr(i);k<s.getRowPtr(i+1);k++) {
					size_t col = s.getCol(k);
					RealType val = s.getValue(k);
					if (col==i) assert(isAlmostZero(val-1.0));
					else assert(isAlmostZero(val));
				}
			}
		}

		bool isAlmostZero(const RealType& x) const
		{
			return (fabs(x)<1e-6);
		}

		void split(SparseMatrixType& matrixA,SparseMatrixType& matrixB,const SparseMatrixType& matrix) const
		{
			size_t counter = 0;
			matrixA.resize(plusSector_);
			for (size_t i=0;i<plusSector_;i++) {
				matrixA.setRow(i,counter);
				for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
					size_t col = matrix.getCol(k);
					RealType val = matrix.getValue(k);
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
					RealType val = matrix.getValue(k);
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

		PsimagLite::ProgressIndicator progress_;
		SparseMatrixType transform_;
		size_t plusSector_;
//		SparseMatrixType s_;
	}; // class ReflectionSymmetry
} // namespace Dmrg

#endif  // REFLECTION_SYMM_H
