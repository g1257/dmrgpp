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

/*! \file ModelCommon.h
 *
 *  Common functions for most models
 *
 */
 
#ifndef MODEL_COMMON_H
#define MODEL_COMMON_H

#include "VerySparseMatrix.h"
#include "IoSimple.h"
#include "Basis.h"

namespace Dmrg {
	//! Common functions for various models
	
	template<typename ModelHelperType,
	typename SparseMatrixType,
 	typename DmrgGeometryType,
	typename LinkProductType,
	typename SharedMemoryType>
	class ModelCommon  {
		public:	
			typedef typename ModelHelperType::RealType RealType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef VerySparseMatrix<SparseElementType> VerySparseMatrixType;
			
			ModelCommon(int DEGREES_OF_FREEDOM,const DmrgGeometryType& dmrgGeometry)
			: dof_(DEGREES_OF_FREEDOM),dmrgGeometry_(dmrgGeometry)
			{
			}

			//! Let H be the hamiltonian of the FeAs model for basis1 and partition m consisting of the external product
			//! of basis2 \otimes basis3
			//! This function does x += H*y
			template<typename SomeVectorType>
			void matrixVectorProduct(
					SomeVectorType& x,
     					const SomeVectorType& y,
	  				ModelHelperType const &modelHelper) const
			{
				//! contribution to Hamiltonian from current system
				modelHelper.hamiltonianLeftProduct(x,y);
				//! contribution to Hamiltonian from current envirnoment
				modelHelper.hamiltonianRightProduct(x,y);
				//! contribution to Hamiltonian from connection system-environment
				hamiltonianConnectionProduct(x,y,modelHelper);
			}

			//! Add Hamiltonian connection between basis2 and basis3 in the orderof basis1 for symmetry block m 
			void addHamiltonianConnection(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
			{
				size_t n = matrix.rank();
				VerySparseMatrixType vsm(n);
				addHamiltonianConnection(vsm,modelHelper);
				SparseMatrixType matrixBlock;
				matrixBlock = vsm;
				matrix += matrixBlock;
			}
			
			//! Add Hamiltonian connection between basis2 and basis3 in the orderof basis1 for symmetry block m 
			void addHamiltonianConnection(VerySparseMatrixType& matrix,
					const ModelHelperType& modelHelper) const
			{
				size_t n=modelHelper.basis1().block().size();
				
				VerySparseMatrixType matrix2(matrix.rank());
				
				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
						SparseMatrixType matrixBlock(matrix.rank(),matrix.rank());
						int res = hamiltonianConnection(i,j,
								&matrixBlock,modelHelper);
						if (res<0) continue;
						VerySparseMatrixType vsm(matrixBlock);
						matrix2+=vsm;
					}
				}
				matrix += matrix2;
			}

			//! Let H_m be the Hamiltonian connection between basis2 and basis3 in the orderof basis1 for block m 
			//! Then this function does x+= H_m *y
			void hamiltonianConnectionProduct(std::vector<SparseElementType> &x,std::vector<SparseElementType> const &y,
				const ModelHelperType& modelHelper) const
			{
				size_t n=modelHelper.basis1().block().size();

				SparseMatrixType matrix;

				typename LinkProductType::LinkProductStructType lps;

				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
						if (i==j) continue; 
						hamiltonianConnection(i,j,0,modelHelper,&lps);
					}
				}
				size_t total = lps.isaved.size();

				LinkProductType pfh(dof_,modelHelper,lps,x,y);
				SharedMemoryType pthreads;
				pthreads.loopCreate(total,pfh);

			}
			
			//! Return H, the hamiltonian of the model for basis1 and partition m consisting of the external product
			//! of basis2 \otimes basis3
			//! Note: Used only for debugging purposes
			void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
			{
				SparseMatrixType matrixBlock;
	
				//! contribution to Hamiltonian from current system
				modelHelper.calcHamiltonianPart(matrixBlock,true);
				matrix = matrixBlock;
				psimag::Matrix<SparseElementType> fm;
				crsMatrixToFullMatrix(fm,matrix);
				std::vector<RealType> e(matrix.rank());
				utils::diag(fm,e,'N');
				std::cerr<<"eSystem[0]="<<e[0]<<"\n";
	
				//! contribution to Hamiltonian from current envirnoment
				modelHelper.calcHamiltonianPart(matrixBlock,false);
				matrix += matrixBlock;
				crsMatrixToFullMatrix(fm,matrixBlock);
				utils::diag(fm,e,'N');
				std::cerr<<"eEnv[0]="<<e[0]<<"\n";
				matrixBlock.clear();
	
				VerySparseMatrixType vsm(matrix);
				addHamiltonianConnection(vsm,modelHelper);
	
				matrix = vsm;
			}


		private:
			size_t dof_;
			const DmrgGeometryType& dmrgGeometry_;
			
			int hamiltonianConnection(size_t i, size_t j,
					SparseMatrixType* matrixBlock,
					const ModelHelperType& modelHelper,
     					typename LinkProductType::LinkProductStructType* lps=0) const
			{
				int flag = -1;
				size_t ind = modelHelper.basis1().block()[i];
				size_t jnd = modelHelper.basis1().block()[j];
				size_t type = dmrgGeometry_.connectionKind(ind,jnd);
				if (type==ProgramGlobals::SYSTEM_SYSTEM || 
					type==ProgramGlobals::ENVIRON_ENVIRON) return flag;
				
				for (size_t term=0;term<dmrgGeometry_.terms();term++) {
					for (size_t dofs=0;dofs<LinkProductType::dofs();dofs++) {
						std::pair<size_t,size_t> edofs = LinkProductType::edofs(dofs,term);
						SparseElementType tmp = dmrgGeometry_(ind,edofs.first,jnd,edofs.second,term);
				
						if (tmp==0.0) continue;
						
						flag += (term  + dofs + 1);
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
									modelHelper,term,dofs);
							*matrixBlock += mBlock;
						}
					}
				}
				return flag;
			}
	};     //class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
