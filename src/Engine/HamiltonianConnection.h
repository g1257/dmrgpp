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
/** \ingroup DMRG */
/*@{*/

/*! \file HamiltonianConnection.h
 *
 * DOC TBW FIXME
 */
#ifndef HAMILTONIAN_CONNECTION_H
#define HAMILTONIAN_CONNECTION_H

#include "Utils.h"
#include "LinkProductStruct.h"

namespace Dmrg {
	
	template<typename GeometryType,typename ModelHelperType,typename LinkProductType>
	class HamiltonianConnection {
		public:
			typedef typename ModelHelperType::RealType RealType;
			typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef LinkProductStruct<SparseElementType> LinkProductStructType;
			typedef typename ModelHelperType::LinkType LinkType;
			typedef std::pair<size_t,size_t> PairType;
			
			HamiltonianConnection(
				const GeometryType& geometry,
				const ModelHelperType& modelHelper,
				const LinkProductStructType* lps = 0,
				std::vector<SparseElementType>* x = 0,
				const std::vector<SparseElementType>* y = 0)
			: lps_(*lps),x_(*x),y_(*y),geometry_(geometry),modelHelper_(modelHelper),
				systemBlock_(modelHelper.basis2().block()),
				envBlock_(modelHelper.basis3().block()),
				smax_(*std::max_element(systemBlock_.begin(),systemBlock_.end())),
				emin_(*std::min_element(envBlock_.begin(),envBlock_.end()))
			{
			}
			
			bool compute(size_t i, size_t j,SparseMatrixType* matrixBlock,
     					LinkProductStructType* lps=0) const
			{
				
				bool flag=false;
				size_t ind = modelHelper_.basis1().block()[i];
				size_t jnd = modelHelper_.basis1().block()[j];
				//throw std::runtime_error("system block is not correct here, think finite algorithm!!!\n"); 
				if (!geometry_.connected(smax_,emin_,ind,jnd)) return flag;
				size_t type = geometry_.connectionKind(smax_,ind,jnd);
				
				if (type==ProgramGlobals::SYSTEM_SYSTEM || 
					type==ProgramGlobals::ENVIRON_ENVIRON) return flag;
				
				for (size_t term=0;term<geometry_.terms();term++) {
					for (size_t dofs=0;dofs<LinkProductType::dofs(term);dofs++) {
						std::pair<size_t,size_t> edofs = LinkProductType::connectorDofs(term,dofs);
						SparseElementType tmp = geometry_(smax_,emin_,
								ind,edofs.first,jnd,edofs.second,term);
				
						if (tmp==0.0) continue;
						
						flag = true;
						// if .. else here is inefficient FIXME
						//std::cerr<<"Adding "<<i<<" "<<j<<" term"<<term<<" dofs="<<dofs<<" value="<<tmp<<"\n";
						if (lps!=0) {
							lps->isaved.push_back(i);
							lps->jsaved.push_back(j);
							//lps->dof1saved.push_back(dof1);
							//lps->dof2saved.push_back(dof2);
							lps->typesaved.push_back(type);
							lps->tmpsaved.push_back(tmp);
							lps->termsaved.push_back(term);
							lps->dofssaved.push_back(dofs);
						} else {
							SparseMatrixType mBlock;
							calcBond(mBlock,i,j,type,tmp,term,dofs);
							*matrixBlock += mBlock;
						}
					}
				}
				return flag;
			}
			
			void thread_function_(size_t threadNum,size_t blockSize,pthread_mutex_t* myMutex)
			{
				std::vector<SparseElementType> xtemp(x_.size(),0);
				//for (size_t i=0;i<xtemp.size();i++) xtemp[i]=0;
				for (size_t p=0;p<blockSize;p++) {
					size_t ix = threadNum * blockSize + p;
					if (ix>=lps_.isaved.size()) break;
					size_t i=lps_.isaved[ix];
					size_t j=lps_.jsaved[ix];
					//size_t dof1=lps_.dof1saved[ix];
					//size_t dof2=lps_.dof2saved[ix];
					size_t type=lps_.typesaved[ix];
					size_t term = lps_.termsaved[ix];
					size_t dofs = lps_.dofssaved[ix];
					SparseElementType tmp=lps_.tmpsaved[ix];
					linkProduct(xtemp,y_,i,j,type,tmp,term,dofs);
					
				}
				if (myMutex) pthread_mutex_lock( myMutex);
				for (size_t i=0;i<x_.size();i++) x_[i]+=xtemp[i];
				if (myMutex) pthread_mutex_unlock( myMutex );
			}

			
		private:
			//! Adds a connector between system and environment
			size_t calcBond(
				SparseMatrixType &matrixBlock,
    				size_t i,
				size_t j,
    				size_t type,
				const SparseElementType& valuec,
				size_t term,
    				size_t dofs) const
			{
				int offset = modelHelper_.basis2().block().size();
				PairType ops;
				std::pair<char,char> mods('N','C');
				size_t fermionOrBoson=ProgramGlobals::FERMION,angularMomentum=0,category=0;
				RealType angularFactor=0;
				bool isSu2 = modelHelper_.isSu2();
				SparseElementType value = valuec;
				LinkProductType::valueModifier(value,term,dofs,isSu2);
				LinkProductType::setLinkData(term,dofs,isSu2,fermionOrBoson,
						ops,mods,angularMomentum,angularFactor,category);
				LinkType link(i,j,type, value,dofs,
					      fermionOrBoson,ops,mods,angularMomentum,angularFactor,category);
				if (link.type==ProgramGlobals::SYSTEM_ENVIRON) {
						
					const SparseMatrixType& A=
						modelHelper_.getReducedOperator(link.mods.first,link.site1,
							link.ops.first,ModelHelperType::System);
					const SparseMatrixType& B=
						modelHelper_.getReducedOperator(link.mods.second,link.site2-offset,
							link.ops.second,ModelHelperType::Environ);
					modelHelper_.fastOpProdInter(A,B,matrixBlock,link);
				} else {
// 						
					if (link.type!=ProgramGlobals::ENVIRON_SYSTEM) std::cerr<<"EEEEEEEEEEEERRRRRRRRRRRRRRROR\n";
					const SparseMatrixType& A=
						modelHelper_.getReducedOperator(link.mods.first,link.site1-offset,
							link.ops.first,ModelHelperType::Environ);
					const SparseMatrixType& B=
						modelHelper_.getReducedOperator(link.mods.second,link.site2,
							link.ops.second,ModelHelperType::System);
					modelHelper_.fastOpProdInter(A,B,matrixBlock,link);
				}
				
				return matrixBlock.nonZero();
				
			}

			//! Computes x+=H_{ij}y where H_{ij} is a Hamiltonian that connects system and environment 
			void linkProduct(std::vector<SparseElementType> &x,std::vector<SparseElementType> const &y,
						size_t i,size_t j,size_t type,
				const SparseElementType  &valuec,size_t term,size_t dofs)  const
			{
				int offset =modelHelper_.basis2().block().size();
				std::pair<size_t,size_t> ops;
				std::pair<char,char> mods('N','C');
				size_t fermionOrBoson=ProgramGlobals::FERMION,angularMomentum=0,category=0;
				RealType angularFactor=0;
				bool isSu2 = modelHelper_.isSu2();
				LinkProductType::setLinkData(term,dofs,isSu2,
						fermionOrBoson,ops,mods,angularMomentum,angularFactor,category);
				SparseElementType value = valuec;
				LinkProductType::valueModifier(value,term,dofs,isSu2);
				LinkType link(i,j,type, value,dofs,
					      fermionOrBoson,ops,mods,angularMomentum,angularFactor,category);
				if (type==ProgramGlobals::SYSTEM_ENVIRON) {
					
					//A=modelHelper.basis2().getOperator(i,sigma);
					const SparseMatrixType& A=modelHelper_.getReducedOperator(mods.first,i,
							link.ops.first,ModelHelperType::System);
					//B=modelHelper.getTcOperator(dof_*(j-offset)+sigma2,ModelHelperType::Environ);
					const SparseMatrixType& B=modelHelper_.getReducedOperator(mods.second,j-offset,
							link.ops.second,ModelHelperType::Environ);
					modelHelper_.fastOpProdInter(x,y,A,B,link);
						
				} else {		
					if (type!=ProgramGlobals::ENVIRON_SYSTEM) std::cerr<<"EEEEEEEEEEEERRRRRRRRRRRRRRROR\n";
						//A=modelHelper.basis3().getOperator(i-offset,sigma);
					const SparseMatrixType& A=modelHelper_.getReducedOperator(mods.first,i-offset,
							link.ops.first,ModelHelperType::Environ);
						
						//B=modelHelper.getTcOperator(j*dof_+sigma2,ModelHelperType::System);
					const SparseMatrixType& B=modelHelper_.getReducedOperator(mods.second,j,
							link.ops.second,ModelHelperType::System);
					modelHelper_.fastOpProdInter(x,y,A,B,link);
						
				}
				
				
			}
			const LinkProductStructType& lps_;
			std::vector<SparseElementType>& x_;
			const std::vector<SparseElementType>& y_;
			const GeometryType& geometry_;
			const ModelHelperType& modelHelper_;
			const typename GeometryType::BlockType& systemBlock_;
			const typename GeometryType::BlockType& envBlock_;
			size_t smax_,emin_;
	}; // class HamiltonianConnection
} // namespace Dmrg 

/*@}*/
#endif // HAMILTONIAN_CONNECTION_H
