/*
Copyright (c) 2009,-2012 UT-Battelle, LLC
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

/*! \file HamiltonianConnection.h
 *
 * DOC TBW FIXME
 */
#ifndef HAMILTONIAN_CONNECTION_H
#define HAMILTONIAN_CONNECTION_H

#include "LinkProductStruct.h"
#include "CrsMatrix.h"

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
			typedef typename GeometryType::AdditionalDataType AdditionalDataType;

			HamiltonianConnection(const GeometryType& geometry,const ModelHelperType& modelHelper,const LinkProductStructType* lps = 0,
			std::vector<SparseElementType>* x = 0,
			const std::vector<SparseElementType>* y = 0)
			: geometry_(geometry),
			  modelHelper_(modelHelper),
			  lps_(*lps),x_(*x),y_(*y),
			  systemBlock_(modelHelper.leftRightSuper().left().block()),
			  envBlock_(modelHelper.leftRightSuper().right().block()),
			  smax_(*std::max_element(systemBlock_.begin(),systemBlock_.end())),
			  emin_(*std::min_element(envBlock_.begin(),envBlock_.end()))
			{}

			bool compute(size_t i,
			             size_t j,
			             SparseMatrixType* matrixBlock,
				     LinkProductStructType* lps,
				     size_t& total) const
			{
				bool flag=false;
				size_t ind = modelHelper_.leftRightSuper().super().block()[i];
				size_t jnd = modelHelper_.leftRightSuper().super().block()[j];

				if (!geometry_.connected(smax_,emin_,ind,jnd)) return flag;
				size_t type = geometry_.connectionKind(smax_,ind,jnd);

				if (type==ProgramGlobals::SYSTEM_SYSTEM || 
					type==ProgramGlobals::ENVIRON_ENVIRON) return flag;

				SparseMatrixType mBlock;

				for (size_t term=0;term<geometry_.terms();term++) {
					geometry_.fillAdditionalData(additionalData_,term,ind,jnd);
					size_t dofsTotal = LinkProductType::dofs(term,additionalData_);
					for (size_t dofs=0;dofs<dofsTotal;dofs++) {
						std::pair<size_t,size_t> edofs = LinkProductType::connectorDofs(term,dofs,additionalData_);
						SparseElementType tmp = geometry_(smax_,emin_,ind,edofs.first,jnd,edofs.second,term);
				
						if (tmp==0.0) continue;
						
						flag = true;
						// if .. else here is inefficient FIXME
						//std::cerr<<"Adding "<<i<<" "<<j<<" term"<<term<<" dofs="<<dofs<<" value="<<tmp<<"\n";
						if (lps!=0) {
							lps->isaved[total]=i;
							lps->jsaved[total]=j;
							//lps->dof1saved.push_back(dof1);
							//lps->dof2saved.push_back(dof2);
							lps->typesaved[total]=type;
							lps->tmpsaved[total]=tmp;
							lps->termsaved[total]=term;
							lps->dofssaved[total]=dofs;
							total++;
						} else {
							calcBond(mBlock,i,j,type,tmp,term,dofs);
							*matrixBlock += mBlock;
						}
					}
				}
				return flag;
			}

			void thread_function_(size_t threadNum,size_t blockSize,size_t total,pthread_mutex_t* myMutex)
			{
				std::vector<SparseElementType> xtemp(x_.size(),0);
				size_t i =0, j = 0, type = 0,term = 0, dofs =0;
				SparseElementType tmp = 0.0;
				for (size_t p=0;p<blockSize;p++) {
					size_t ix = threadNum * blockSize + p;
					if (ix>=total) break;
					prepare(ix,i,j,type,tmp,term,dofs);

					linkProduct(xtemp,y_,i,j,type,tmp,term,dofs);
					
				}
				if (myMutex) pthread_mutex_lock( myMutex);
				for (size_t i=0;i<x_.size();i++) x_[i]+=xtemp[i];
				if (myMutex) pthread_mutex_unlock( myMutex );
			}

			void prepare(size_t ix,size_t& i,size_t& j,size_t& type,SparseElementType& tmp,size_t& term,size_t& dofs) const
			{
				i=lps_.isaved[ix];
				j=lps_.jsaved[ix];
				type=lps_.typesaved[ix];
				term = lps_.termsaved[ix];
				dofs = lps_.dofssaved[ix];
				tmp=lps_.tmpsaved[ix];
				size_t ind = modelHelper_.leftRightSuper().super().block()[i];
				size_t jnd = modelHelper_.leftRightSuper().super().block()[j];
				geometry_.fillAdditionalData(additionalData_,term,ind,jnd);
			}

			LinkType getKron(const SparseMatrixType** A,
					 const SparseMatrixType** B,
					 size_t i,
					 size_t j,
					 size_t type,
					 const SparseElementType& valuec,
					 size_t term,
					 size_t dofs) const
			{
				int offset = modelHelper_.leftRightSuper().left().block().size();
				PairType ops;
				std::pair<char,char> mods('N','C');
				size_t fermionOrBoson=ProgramGlobals::FERMION,angularMomentum=0,category=0;
				RealType angularFactor=0;
				bool isSu2 = modelHelper_.isSu2();
				SparseElementType value = valuec;
				LinkProductType::valueModifier(value,term,dofs,isSu2,additionalData_);
				LinkProductType::setLinkData(term,dofs,isSu2,fermionOrBoson,ops,mods,angularMomentum,angularFactor,category,additionalData_);
				LinkType link2(i,j,type, value,dofs,fermionOrBoson,ops,mods,angularMomentum,angularFactor,category);
				size_t sysOrEnv = (link2.type==ProgramGlobals::SYSTEM_ENVIRON) ? ModelHelperType::System : ModelHelperType::Environ;
				size_t envOrSys = (link2.type==ProgramGlobals::SYSTEM_ENVIRON) ? ModelHelperType::Environ : ModelHelperType::System;
				size_t site1Corrected =(link2.type==ProgramGlobals::SYSTEM_ENVIRON) ? link2.site1 : link2.site1-offset;
				size_t site2Corrected =(link2.type==ProgramGlobals::SYSTEM_ENVIRON) ? link2.site2-offset : link2.site2;

				*A = &modelHelper_.getReducedOperator(link2.mods.first,site1Corrected,link2.ops.first,sysOrEnv);
				*B = &modelHelper_.getReducedOperator(link2.mods.second,site2Corrected,link2.ops.second,envOrSys);
//				printFullMatrix(**A,"A");
//				printFullMatrix(**B,"B");
//				std::cout<<"link2.value="<<link2.value<<"i="<<i<<" j="<<j<<" type="<<type<<" valuec="<<valuec<<" term="<<term<<" dofs="<<dofs<<"\n";
				return link2;
			}


			template<typename SomeConcurrencyType,typename SomeOtherConcurrencyType>
			void sync(SomeConcurrencyType& conc,SomeOtherConcurrencyType& conc2)
			{
				conc.reduce(x_,conc2);
			}

		private:

			//! Adds a connector between system and environment
			size_t calcBond(SparseMatrixType &matrixBlock,
					size_t i,
			                size_t j,
					size_t type,
			                const SparseElementType& valuec,
			                size_t term,
					size_t dofs) const
			{
				SparseMatrixType const* A = 0;
				SparseMatrixType const* B = 0;
				LinkType link2 = getKron(&A,&B,i,j,type,valuec,term,dofs);
				modelHelper_.fastOpProdInter(*A,*B,matrixBlock,link2);

				return matrixBlock.nonZero();
			}

			//! Computes x+=H_{ij}y where H_{ij} is a Hamiltonian that connects system and environment 
			void linkProduct(std::vector<SparseElementType> &x,
					 std::vector<SparseElementType> const &y,
					 size_t i,
					 size_t j,
					 size_t type,
					 const SparseElementType &valuec,
					 size_t term,
					 size_t dofs) const
			{
				SparseMatrixType const* A = 0;
				SparseMatrixType const* B = 0;
				LinkType link2 = getKron(&A,&B,i,j,type,valuec,term,dofs);
				modelHelper_.fastOpProdInter(x,y,*A,*B,link2);
			}

			const GeometryType& geometry_;
			const ModelHelperType& modelHelper_;
			const LinkProductStructType& lps_;
			std::vector<SparseElementType>& x_;
			const std::vector<SparseElementType>& y_;
			const typename GeometryType::BlockType& systemBlock_;
			const typename GeometryType::BlockType& envBlock_;
			size_t smax_,emin_;
			mutable AdditionalDataType additionalData_;
	}; // class HamiltonianConnection
} // namespace Dmrg 

/*@}*/
#endif // HAMILTONIAN_CONNECTION_H
