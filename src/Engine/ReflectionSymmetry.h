// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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
#ifndef REFLECTION_SYMMETRY_H
#define REFLECTION_SYMMETRY_H

#include "DmrgBasisWithOperators.h"

/** \ingroup DMRG */
/*@{*/

/*! \file ReflectionSymmetry.h
 *
 *  A class to implement the reflection symmetry present in the inf. algorithm when system and environm are symmetry wrt. reflection
 *
 */

namespace Dmrg { 	
	class ReflectionSymmetry {
	public:
	
		ReflectionSymmetry(bool useReflection) :
			useReflection_(useReflection),
			sqrt2_(sqrt(2.0)),
			reflectionSelves_(0),
			reflectionSector_(911)
		{
		
			createReflectionPermutation(); 
		}
		
		int size(int tmp) const
		{
			
			if (!useReflection_) return tmp;
			if (reflectionSector_==0) return size_t((tmp+reflectionSelves_)/2);
			if (reflectionSector_==1) return size_t((tmp-reflectionSelves_)/2);
			throw std::runtime_error("size() : reflection sector undefined\n");
			return -1; // should never reach here
		}
		
		void getReflectedEigs(
				MatrixElementType& energyTmp,std::vector<MatrixElementType>& tmpVec,
				MatrixElementType energyTmp1,const std::vector<MatrixElementType>& tmpVec1,
				MatrixElementType energyTmp2,const std::vector<MatrixElementType>& tmpVec2) const
		{
			// first set the energy
			energyTmp=energyTmp1;
			int type=0;
			if (energyTmp1>energyTmp2) {
				energyTmp=energyTmp2;
				type=1;
			}
			// now set the eigenvector
			int fullsize= basis1_.partition(m_+1)-basis1_.partition(m_);
			tmpVec.resize(tmpVec1.size()+tmpVec2.size());
			if (tmpVec.size()!=size_t(fullsize)) {
				std::cerr<<"fullsize="<<fullsize<<" but tmpvec.size="<<tmpVec.size()<<"\n"; 
				throw std::runtime_error("getReflectedEigs(): tmpVec.size()!=size_t(fullsize)\n");
			}
			size_t reflectedMinusSize=size_t((reflectionPermutation_.size()-reflectionSelves_)/2);
			if (tmpVec2.size()!=reflectedMinusSize) throw std::runtime_error("getReflectedEigs(): tmpvec2.size!=reflectedMinusSize\n");
			
			for (size_t i=0;i<tmpVec.size();i++) {
				size_t ii=i;
				tmpVec[i]=0;
				MatrixElementType factor = 0.7071067811865474761;
				MatrixElementType sign = 1;
				if (i>=tmpVec1.size()) {
					ii -= tmpVec1.size();
					if (type==1) sign= -1;
				}
				if (type==0) {
					if (i>=tmpVec2.size() && i<tmpVec1.size()) factor =1;
					tmpVec[i]=factor * tmpVec1[ii];
				} else {
					if (i>=tmpVec2.size() && i<tmpVec1.size()) continue;
					tmpVec[i]=factor * tmpVec2[ii] * sign;
				}
			}
			
		}
		
		void setReflectionSymmetry(size_t reflectionSector)
		{
			reflectionSector_=reflectionSector;
			std::cerr<<"Set reflectionSector="<<reflectionSector_;
			//std::cerr<<" size="<<size()<<" m="<<m_<<"\n";
			//for (int i=0;i<basis1_.partition(m_+1)-basis1_.partition(m_);i++) 
			//	std::cerr<<"qn["<<i<<"]="<<basis1_.qn(i+basis1_.partition(m_))<<" reflx="<<getReflectedState(i)<<"\n";	
		}
		
		void printFullMatrix(const SparseMatrixType& matrix) const
                {
                        size_t n = matrix.getSize();
			SparseMatrixType matrix2=matrix;
			matrix2=reflectionPermute(matrix);
                        std::cerr<<"Matrix for lanczos of size n="<<n<<"\n";
                        for (size_t r=0;r<n;r++) {
                                for (size_t j=0;j<n;j++) {
                                        std::cerr<<matrix2(r,j)<<"\t";
                                }
                                std::cerr<<"\n";
                        }
                }	
		
		// reflection functions below
		bool outsideReflectionBounds(int i) const
		{
			if (!useReflection_) return false;
			int reflectedMinusSize=size_t((reflectionPermutation_.size()-reflectionSelves_)/2);
			int reflectedPlusSize=reflectionPermutation_.size()-reflectedMinusSize;
			int ii = reflectionPermutation_[i]-reflectionSector_*reflectedPlusSize;
			if (ii>=size() || ii<0) return true;
			return false;
		}
		
		void elementMultiplication(MatrixElementType value, std::vector<MatrixElementType>& x,
				const std::vector<MatrixElementType>& y, int i,int j) const
		{
			if (!useReflection_) {
				x[i]+= value*y[j];
				return;
			}

			MatrixElementType tmp=1;
			int reflectedMinusSize=size_t((reflectionPermutation_.size()-reflectionSelves_)/2);
			int reflectedPlusSize=reflectionPermutation_.size()-reflectedMinusSize;
			int jj = reflectionPermutation_[j]-reflectionSector_*reflectedPlusSize;
			int ii = reflectionPermutation_[i]-reflectionSector_*reflectedPlusSize;
			if (jj<0) {
				jj += size();
				if (reflectionSector_==1) tmp = -1;
				if (jj>=reflectedMinusSize || ii>=reflectedMinusSize) return;
			}
			//if (reflectionSector_==1 && reflectedMinusSize!=size_t(size())) throw std::runtime_error("internal at elementMultiplication\n");
			if (jj>=size()) {
				jj -= size();
				if (reflectionSector_==1) tmp = -1;
				if (jj>=reflectedMinusSize || ii>=reflectedMinusSize) return;
			}
			
			if (reflectionSector_==0) {
				if (jj>=reflectedMinusSize && ii<reflectedMinusSize) tmp *= sqrt2_;
				if (ii>=reflectedMinusSize && jj<reflectedMinusSize) tmp *= sqrt2_;
			}
			//if (reflectionSelves_==0 && tmp==sqrt2_) throw std::runtime_error("internal at elementMultiplication\n");
			x[ii] += value*y[jj]*tmp;
		}
		
	private:
		MatrixElementType sqrt2_;
		size_t reflectionSelves_;
		size_t reflectionSector_;
		std::vector<size_t> reflectionPermutation_;
		
		SparseMatrixType reflectionPermute(const SparseMatrixType& matrix) const
		{
			int n = matrix.getSize();	
			psimag::Matrix<MatrixElementType> matrix2(n,n);
			for (int i=0;i<n;i++) {
				for (int j=0;j<n;j++) {
					matrix2(reflectionPermutation_[i],reflectionPermutation_[j])=matrix(i,j);
				}
			}
			SparseMatrixType matrix3(matrix2);
			return matrix3;
		}	
		
		void createReflectionPermutation()
		{
			
			reflectionPermutation_.resize(basis1_.partition(m_+1)-basis1_.partition(m_));
			
			if (!useReflection_) {
				for (size_t i=0;i<reflectionPermutation_.size();i++) reflectionPermutation_[i]=i;
				return;
			}
				
			if (basis2_.size()!=basis3_.size()) throw std::runtime_error("DmrgModelHelper::createReflectionPermutation: not reflection symmetry\n");
			
			//checkReflectionFunction();
	
			size_t j=0;	
			std::vector<size_t> rvector,selves;
			for (size_t i=0;i<reflectionPermutation_.size();i++) {
				
				size_t k = getReflectedState(i);
				if (utils::isInVector(rvector,i)>=0) continue;
				if (utils::isInVector(rvector,k)>=0) continue;
				if (k==i) {
					std::cerr<<"FFFFFFFFFFFFFF"<<k<<"\n";
					if (utils::isInVector(selves,i)<0) selves.push_back(i);
					continue;
				} 
				rvector.push_back(k);
				
				reflectionPermutation_[j]=i;
				j++;
			}
			reflectionSelves_=selves.size();
			for (size_t i=0;i<selves.size();i++) {
				std::cerr<<"SELVES["<<i<<"]="<<selves[i]<<"\n";
				reflectionPermutation_[j]=selves[i];
				j++;
			}
			for (size_t i=0;i<rvector.size();i++) {
				reflectionPermutation_[j]=rvector[i];
				j++;
			}
			invert(reflectionPermutation_);	
			//std::cerr<<"Resetting reflection permutation here!!!!!!!\n";	
			//for (size_t i=0;i<reflectionPermutation_.size();i++) reflectionPermutation_[i]=i;	
			utils::vectorPrint(reflectionPermutation_,"reflectionPermutation_",std::cerr);		
			std::cerr<<"Found selves="<<reflectionSelves_<<"\n";
				
		}
		
		void checkReflectionFunction()
		{
			for (size_t i=0;i<reflectionPermutation_.size();i++) {
				size_t j = getReflectedState(i);
				if (getReflectedState(j)!=i) {
					throw std::runtime_error("Reflection^2 is not the identity!!\n");
				}
			}
		}		
		
		template<typename T>		
		void invert(std::vector<T>& v)
		{
			std::vector<T> w=v;
			for (size_t i=0;i<v.size();i++) v[w[i]]=i;
		}
	
		// get the state reflected 
		size_t getReflectedState(size_t i)
		{
			int smallHilbertSpaceSize=4; 
			int nsPrev = int(basis2_.size()/smallHilbertSpaceSize);
			int nePrev = smallHilbertSpaceSize;
			if (basis2_.size() % smallHilbertSpaceSize !=0) {
				throw std::runtime_error("basis2_.size() not divisible by smallHilbertSpaceSize\n");
			}
			if (basis3_.size() % smallHilbertSpaceSize !=0) {
                                throw std::runtime_error("basis3.size() not divisible by smallHilbertSpaceSize\n");
                        }

			int alpha0,alpha1, beta0,beta1;
			utils::getCoordinates(alpha0,beta0,basis2_.permutation(alpha_[i]),nsPrev);
			utils::getCoordinates(alpha1,beta1,basis3_.permutation(beta_[i]),nePrev);
			int alphaNew=basis2_.permutationInverse(beta1+alpha1*nsPrev);
			int betaNew= basis3_.permutationInverse(beta0+alpha0*nePrev);
			
			int tmp= buffer_[alphaNew][betaNew];
			if (tmp<0) throw std::runtime_error("DmrgModelHelper::getReflectedState(): reflection of state not in same symmetry block\n");
			return tmp;
		}

		size_t internalRefl(size_t i,size_t n)
		{
			size_t y = size_t(i/n);
			size_t x = i % n; 
			return y + x*n;
		}
	};
} // namespace Dmrg
/*@}*/

#endif

