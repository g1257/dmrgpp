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

/*! \file HamiltonianConnection.h
 *
 * DOC TBW FIXME
 */
#ifndef HAMILTONIAN_CONNECTION_H
#define HAMILTONIAN_CONNECTION_H

#include "Utils.h"

namespace Dmrg {
	
	template<typename GeometryType,typename ModelHelperType,typename LinkProductType>
	class HamiltonianConnection {
		public:
			typedef typename ModelHelperType::RealType RealType;
			typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef typename LinkProductType::LinkProductStructType LinkProductStructType;
			
			HamiltonianConnection(const GeometryType& geometry,const ModelHelperType& modelHelper)
			: geometry_(geometry),modelHelper_(modelHelper),
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
					for (size_t dofs=0;dofs<LinkProductType::dofs();dofs++) {
						std::pair<size_t,size_t> edofs = LinkProductType::edofs(dofs,term);
						SparseElementType tmp = geometry_(smax_,emin_,
								ind,edofs.first,jnd,edofs.second,term);
				
						if (tmp==0.0) continue;
						
						flag = true;
						// if .. else here is inefficient FIXME
						//std::cerr<<"Adding "<<i<<" "<<j<<" "<<connectionType<<" value="<<tmp<<"\n";
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
							LinkProductType::calcBond(i,j,type,tmp,mBlock,
									modelHelper_,term,dofs);
							*matrixBlock += mBlock;
						}
					}
				}
				return flag;
			}
			
			
		private:
			const GeometryType& geometry_;
			const ModelHelperType& modelHelper_;
			const typename GeometryType::BlockType& systemBlock_;
			const typename GeometryType::BlockType& envBlock_;
			size_t smax_,emin_;
	}; // class HamiltonianConnection
} // namespace Dmrg 

/*@}*/
#endif // HAMILTONIAN_CONNECTION_H
