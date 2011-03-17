// BEGIN LICENSE BLOCK
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
#ifndef RIGHT_LEFT_LOCAL_H
#define RIGHT_LEFT_LOCAL_H

#include "BLAS.h"

/** \ingroup DMRG */
/*@{*/

/*! \file RightLeftLocal.h
 *
 *  A class to contain state information about the Hamiltonian to help with the calculation of x+=Hy
 *
 */

namespace Dmrg { 	
	template<typename BasisType,typename BasisWithOperatorsType,typename SparseMatrixType>
	class RightLeftLocal {
	public:	
		typedef typename SparseMatrixType::value_type MatrixElementType;
		typedef psimag::Matrix<MatrixElementType> MatrixType;
		
		RightLeftLocal(int m,const BasisType& basis1,  const BasisWithOperatorsType& basis2,
			 const BasisWithOperatorsType& basis3,size_t orbitals,bool useReflection=false) 
		:
			m_(m),
			basis1_(basis1),
			basis2_(basis2),
			basis3_(basis3),
			alpha_(basis1_.size()),
			beta_(basis1_.size()),
			leftPermInv_(basis2.size()),
			rightPermInv_(basis3.size())
		{
			init();
			createAlphaAndBeta();
		}

		~RightLeftLocal()
		{
			//for (size_t i=0;i<bMatrix_.size();i++) delete bMatrix_[i];
			//for (size_t i=0;i<aMatrix_.size();i++) delete aMatrix_[i];
		}

		//! Does x+= (AB)y, where A belongs to pSprime and B  belongs to pEprime or viceversa (inter)
		//! Has been changed to accomodate for reflection symmetry
		void fastOpProdInter(	std::vector<MatrixElementType>  &x,
					std::vector<MatrixElementType>  const &y,
					SparseMatrixType const &A,
					SparseMatrixType const &B,
					int type,
					MatrixElementType  &hop,
					bool operatorsAreFermions=true,size_t angularMomentum=1,MatrixElementType angularSign= -1.0,size_t category=0,bool dummy2=false) const
		{
			int const SystemEnviron=1,EnvironSystem=2;
			int fermionSign =  (operatorsAreFermions) ? -1 : 1;
			
			
			if (type==EnvironSystem)  {
				MatrixElementType hop2 =hop*fermionSign;
				fastOpProdInter(x,y,B,A,SystemEnviron,hop2,operatorsAreFermions);
				return;
			}
			size_t leftSize = leftPerm_.size();
			size_t rightSize = rightPerm_.size();
			//static const std::vector<MatrixElementType>* yAddress = 0;
			
			//if (yAddress!=&y) {
				preparePhi(yMatrix_,y);
				prepareB(bMatrix_,B);
				prepareA(aMatrix_,A,operatorsAreFermions);
			//	yAddress = &y;
			//}
			
			MatrixType* bm = &bMatrix_;
			/*int ib = PsimagLite\:\:isInVector(addressesB_,&B);
			
			if (ib<0) {
				bm = new MatrixType(rightSize,rightSize);
				prepareB(*bm,B);
				bMatrix_.push_back(bm);
				addressesB_.push_back(&B);
			} else {
				bm =  bMatrix_[ib];
			}
			
			int ia = PsimagLite\:\:isInVector(addressesA_,&A);*/
			MatrixType* am = &aMatrix_;
			/*if (ia<0) {
				am = new MatrixType(leftSize,leftSize);
				prepareA(*am,A,operatorsAreFermions);
				aMatrix_.push_back(am);
				addressesA_.push_back(&A);
			} else {
				am =  aMatrix_[ia];
			}*/
			
			
			//! multiply all here:
			
			psimag::BLAS::GEMM('N','C',rightSize,leftSize,rightSize,hop,
					   &(bm->operator()(0,0)),rightSize,&(yMatrix_(0,0)),leftSize,0.0,&(cMatrix_(0,0)),rightSize);
			psimag::BLAS::GEMM('N','C',leftSize,rightSize,leftSize,1.0,
					   &(am->operator()(0,0)),leftSize,&(cMatrix_(0,0)),rightSize,0.0,&(tmpMatrix_(0,0)),leftSize);
			
			//! revert order
			unpreparePhi(x,tmpMatrix_);
		}

	private:
		int m_;
		const BasisType&  basis1_;
		const BasisWithOperatorsType& basis2_;
		const BasisWithOperatorsType& basis3_;
		std::vector<size_t> alpha_,beta_;
		std::vector<size_t> leftPermInv_,rightPermInv_;
		std::vector<size_t> leftPerm_,rightPerm_;
		mutable MatrixType bMatrix_;
		mutable MatrixType aMatrix_;
		mutable MatrixType cMatrix_,tmpMatrix_,yMatrix_;
		//mutable std::vector<const SparseMatrixType*> addressesA_;
		//mutable std::vector<const SparseMatrixType*> addressesB_;
		
		
		void init()
		{
			size_t ns=basis2_.size();
			size_t ne=basis3_.size();
			int offset = basis1_.partition(m_);
			int total = basis1_.partition(m_+1) - offset;

			for (size_t alphaPrime=0;alphaPrime<ns;alphaPrime++) {
				for (size_t betaPrime=0;betaPrime<ne;betaPrime++) {	
					int tmp =basis1_.permutationInverse(alphaPrime + betaPrime*ns) - offset;
					if (tmp>=total || tmp<0) continue;
					int x = PsimagLite\:\:isInVector(leftPerm_,alphaPrime);
					if (x<0) leftPerm_.push_back(alphaPrime);
					int y = PsimagLite\:\:isInVector(rightPerm_,betaPrime);
					if (y<0) rightPerm_.push_back(betaPrime);
				}
			}
			
			
			for (size_t i=0;i<rightPerm_.size();i++) rightPermInv_[rightPerm_[i]]=i;
			
			for (size_t i=0;i<leftPerm_.size();i++) leftPermInv_[leftPerm_[i]]=i;
			
			size_t leftSize = leftPerm_.size();
			size_t rightSize = rightPerm_.size();
			
			yMatrix_.resize(leftSize,rightSize);
			cMatrix_.resize(rightSize,leftSize);
			tmpMatrix_.resize(leftSize,rightSize);
			aMatrix_.resize(leftSize,leftSize);
			bMatrix_.resize(rightSize,rightSize);
			
		}
		
		void preparePhi(MatrixType& m,std::vector<MatrixElementType>  const &v) const
		{
			int offset = basis1_.partition(m_);
			int total = basis1_.partition(m_+1) - offset;
			/*for (size_t i=0;i<leftPerm_.size();i++) {
				size_t x = leftPerm_[i];
				for (size_t j=0;j<rightPerm_.size();j++) {
					size_t y = rightPerm_[j];
					int ii = basis1_.permutationInverse(x+y*basis2_.size())-offset;
					if (ii<0 || ii>=total) continue;
					m(i,j) = v[ii];
				}
			}*/
			
			//size_t ns = basis2_.size();
			for (int i=0;i<total;i++) {
				//size_t alpha,beta;
				//utils::getCoordinates(alpha,beta,basis1_.permutation(i+offset),ns);
				m(leftPermInv_[alpha_[i]],rightPermInv_[beta_[i]])=v[i];
			}
		}
		
		void unpreparePhi(std::vector<MatrixElementType>& v,MatrixType& m) const
		{
			int offset = basis1_.partition(m_);
			int total = basis1_.partition(m_+1) - offset;
			/*for (size_t i=0;i<leftPerm_.size();i++) {
				size_t x = leftPerm_[i];
				for (size_t j=0;j<rightPerm_.size();j++) {
					size_t y = rightPerm_[j];
					int ii = basis1_.permutationInverse(x+y*basis2_.size())-offset;
					if (ii<0 || ii>=total) continue;
					//MatrixElementType a =v[ii];
					//a+=2.0;
					 v[ii] += m(i,j);
				}
			}*/
			
			size_t ns = basis2_.size();
			for (int i=0;i<total;i++) {
				size_t alpha,beta;
				utils::getCoordinates(alpha,beta,basis1_.permutation(i+offset),ns);
				v[i] =m(leftPermInv_[alpha_[i]],rightPermInv_[beta_[i]]);
			}
		}
		
		void prepareB(MatrixType& m,SparseMatrixType const &B) const
		{
			for (size_t i=0;i<rightPerm_.size();i++) {
				size_t x = rightPerm_[i];
				for (size_t j=0;j<rightPerm_.size();j++) {
					size_t y = rightPerm_[j];
					m(i,j) = B(x,y);
				}
			}
		}
		
		void prepareA(MatrixType& m,SparseMatrixType const &A,
			      bool operatorsAreFermions) const
		{
			int fermionSign =  (operatorsAreFermions) ? -1 : 1;
			for (size_t i=0;i<leftPerm_.size();i++) {
				size_t x = leftPerm_[i];
				MatrixElementType tmp = basis2_.fermionicSign(x,fermionSign);
				for (size_t j=0;j<leftPerm_.size();j++) {
					size_t y = leftPerm_[j];
					m(i,j) = A(x,y)*tmp;
				}
			}
		}
		
		void createAlphaAndBeta()
		{
			size_t ns=basis2_.size();
			int offset = basis1_.partition(m_);
			int total = basis1_.partition(m_+1) - offset;

			for (int i=0;i<total;i++) {
				// row i of the ordered product basis
				utils::getCoordinates(alpha_[i],beta_[i],basis1_.permutation(i+offset),ns);
			}
		}
		
	}; // class RightLeftLocal
} // namespace Dmrg
/*@}*/

#endif

