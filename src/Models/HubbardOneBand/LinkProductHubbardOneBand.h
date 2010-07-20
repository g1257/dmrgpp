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
/** \ingroup DMRG */
/*@{*/

/*! \file LinkProductHubbardOneBand.h
 *
 *  A class to represent a generic solver for the Dmrg method
 *
 */
#ifndef LINK_PRODUCT_HubbardOneBand_H
#define LINK_PRODUCT_HubbardOneBand_H

#include "LinkProductStruct.h"

namespace Dmrg {
	
	
	
	template<typename ModelHelperType>
	class LinkProductHubbardOneBand {
			typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
			
		public:
			typedef typename ModelHelperType::RealType RealType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef LinkProductStruct<SparseElementType> LinkProductStructType;
			
			LinkProductHubbardOneBand(size_t dof,
				   const ModelHelperType& modelHelper,
				   const LinkProductStructType& lps,
				    std::vector<SparseElementType>& x,
				   const std::vector<SparseElementType>& y) : dof_(dof),modelHelper_(modelHelper),lps_(lps),
					x_(x),y_(y)
			{
				
			}
			
			void thread_function_(size_t threadNum,size_t blockSize,pthread_mutex_t* myMutex)
			{
				std::vector<SparseElementType> xtemp(x_.size(),0);
				//for (size_t i=0;i<xtemp.size();i++) xtemp[i]=0;
				for (size_t p=0;p<blockSize;p++) {
					size_t ix = threadNum * blockSize + p;
					size_t i=lps_.isaved[ix];
					size_t j=lps_.jsaved[ix];
					//size_t dof1=lps_.dof1saved[ix];
					//size_t dof2=lps_.dof2saved[ix];
					int type=lps_.typesaved[ix];
					size_t term = lps_.termsaved[ix];
					size_t dofs = lps_.dofssaved[ix];
					SparseElementType tmp=lps_.tmpsaved[ix];
					linkProduct(xtemp,y_,i,j,type,tmp,modelHelper_,term,dofs);
					
				}	
				if (myMutex) pthread_mutex_lock( myMutex);
				for (size_t i=0;i<x_.size();i++) x_[i]+=xtemp[i];
				if (myMutex) pthread_mutex_unlock( myMutex );
				
			} 
			
			//! Adds a tight-binding bond between system and environment
			static size_t calcBond(size_t i,size_t j,size_t type,
				SparseElementType  &val,
				SparseMatrixType &matrixBlock,
				const ModelHelperType& modelHelper,size_t term,size_t dofs)
			{
				size_t sigma = dofs;
				size_t sigma2 = dofs;
				//if (sigma!=sigma2) return 0;
				//int const SystemEnviron=1,EnvironSystem=2;
				SparseMatrixType A,B;
				//int k=sigma + i*pSprime.dof();
				//int k2=sigma2 + j*pSprime.dof();
				int offset = modelHelper.basis2().block().size();
				size_t angularMomentum=1;
				RealType angularFactor  = 1.0;
				if (sigma==1) angularFactor  = -1.0;
				
				if (type==ProgramGlobals::SYSTEM_ENVIRON) {
					
// 					A=modelHelper.basis2().getOperator(i,sigma).data;
// 					transposeConjugate(B,modelHelper.basis3().getOperator(j-offset,sigma).data);
// 					angularMomentum=0;
//					modelHelper.fastOpProdInter(A,B,type,val,matrixBlock,true);
					const SparseMatrixType& A=modelHelper.getReducedOperator('N',i,sigma,ModelHelperType::System);
					const SparseMatrixType& B=modelHelper.getReducedOperator('C',j-offset,sigma2,ModelHelperType::Environ);
					modelHelper.fastOpProdInter(A,B,type,val,matrixBlock,true,angularMomentum,angularFactor,sigma);
							
					
				} else {
					if (type!=ProgramGlobals::ENVIRON_SYSTEM) std::cerr<<"ERRRRRRRRRRORRRRRRRR\n";
// 					A=modelHelper.basis3().getOperator(i-offset,sigma).data;
// 					transposeConjugate(B,modelHelper.basis2().getOperator(j,sigma).data);
// 					modelHelper.fastOpProdInter(A,B,type,val,matrixBlock,true);
					const SparseMatrixType& A=modelHelper.getReducedOperator('N',i-offset,sigma,ModelHelperType::Environ);
					const SparseMatrixType& B=modelHelper.getReducedOperator('C',j,sigma2,ModelHelperType::System);
					modelHelper.fastOpProdInter(A,B,type,val,matrixBlock,true,angularMomentum,angularFactor,sigma);
					
				}
				
				
				return matrixBlock.nonZero();
				
			}
			
			// used only for testing 
			RealType norm() const
			{
				return std::norm(x_);
			}
			
			static size_t dofs() { return 2; }
			
			// up up and down down are the only connections possible for this model
			static std::pair<size_t,size_t> edofs(size_t dofs,size_t term)
			{
				return std::pair<size_t,size_t>(0,0); // no orbital and no dependence on spin
			}
			
		private:
			size_t dof_;
			const ModelHelperType& modelHelper_;
			const LinkProductStructType& lps_;
			std::vector<SparseElementType>& x_;
			const std::vector<SparseElementType>& y_;
		
			//! Computes x+=H_{ij}y where H_{ij} is a Hamiltonian that connects system and environment 
			void linkProduct(std::vector<SparseElementType> &x,std::vector<SparseElementType> const &y,
					int i,int j,int type,
			     SparseElementType  &val,
			     ModelHelperType const &modelHelper,size_t term,size_t dofs)  const
			{
				size_t sigma = dofs;
				size_t sigma2 = dofs;
				//int const SystemEnviron=1,EnvironSystem=2;
				
				int offset = modelHelper.basis2().block().size();
				size_t angularMomentum=1;
				RealType angularFactor  = 1.0;
				if (sigma==1) angularFactor  = -1.0;
				
				if (type==ProgramGlobals::SYSTEM_ENVIRON) {
					
					//A=modelHelper.basis2().getOperator(i,sigma);
					//transposeConjugate(B,modelHelper.basis3().getOperator(j-offset,sigma));
					const SparseMatrixType& A=modelHelper.getReducedOperator('N',i,sigma,ModelHelperType::System);
					const SparseMatrixType& B=modelHelper.getReducedOperator('C',j-offset,sigma2,ModelHelperType::Environ);
					modelHelper.fastOpProdInter(x,y,A,B,type,val,true,angularMomentum,angularFactor,sigma);
							
				} else {
					if (type!=ProgramGlobals::ENVIRON_SYSTEM) std::cerr<<"ERRRRRRRRRRORRRRRRRR\n";
						//A=modelHelper.basis3().getOperator(i-offset,sigma);
						//transposeConjugate(B,modelHelper.basis2().getOperator(j,sigma));
					const SparseMatrixType& A=modelHelper.getReducedOperator('N',i-offset,sigma,ModelHelperType::Environ);
					const SparseMatrixType& B=modelHelper.getReducedOperator('C',j,sigma2,ModelHelperType::System);
					modelHelper.fastOpProdInter(x,y,A,B,type,val,true,angularMomentum,angularFactor,sigma);
				}
	
				
			}
			
	}; // class LinkProductHubbardOneBand
} // namespace Dmrg
/*@}*/
#endif
