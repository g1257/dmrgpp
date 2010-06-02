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
				int smax,emin;
				
				dmrgGeometry_.findExtremes(smax,emin,modelHelper.basis1().block());
				
				
				VerySparseMatrixType matrix2(matrix.rank());
				
				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
						//if (modelHelper.basis1().block()[i]>modelHelper.basis1().block()[j]) continue; // do not count twice
						for (size_t dof1=0;dof1<dof_;dof1++) {
							for (size_t dof2=0;dof2<dof_;dof2++) {
								//if (i==j) continue; 
								SparseMatrixType matrixBlock(matrix.rank(),matrix.rank());
								int res = hamiltonianConnection(i,dof1,j,dof2,smax,emin,
										&matrixBlock,modelHelper);
								if (res<0) continue;
								VerySparseMatrixType vsm(matrixBlock);
								//matrixBlock.clear();
								matrix2+=vsm;
							}
						}
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
				int smax,emin;

				SparseMatrixType matrix;

				dmrgGeometry_.findExtremes(smax,emin,modelHelper.basis1().block());

				typename LinkProductType::LinkProductStructType lps;

				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
						for (size_t dof1=0;dof1<dof_;dof1++) {
							for (size_t dof2=0;dof2<dof_;dof2++) { 	
								if (i==j) continue; 
								hamiltonianConnection(i,dof1,j,dof2,smax,emin,0,modelHelper,&lps);
								
							}
						}
					}
				}
				size_t total = lps.isaved.size();

				LinkProductType pfh(dof_,modelHelper,lps,x,y);
				SharedMemoryType pthreads;
				pthreads.loopCreate(total,pfh);

			}
			
			//! Return H, the hamiltonian of the FeAs model for basis1 and partition m consisting of the external product
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
				
	
				//! contribution to Hamiltonian from current envirnoment
				modelHelper.calcHamiltonianPart(matrixBlock,false);
				matrix += matrixBlock;
				matrixBlock.clear();
	
				VerySparseMatrixType vsm(matrix);
				addHamiltonianConnection(vsm,modelHelper);
	
				matrix = vsm;
			}


		private:
			size_t dof_;
			const DmrgGeometryType& dmrgGeometry_;
			
			int hamiltonianConnection(size_t i,size_t dof1, size_t j,size_t dof2,
					size_t smax,
     					size_t emin,
					SparseMatrixType* matrixBlock,
					const ModelHelperType& modelHelper,
     					typename LinkProductType::LinkProductStructType* lps=0) const
			{
				int flag = -1;
				for (size_t connectionType=0;connectionType<dmrgGeometry_.connectorValues();connectionType++) {
					int type = dmrgGeometry_.calcConnectorType(modelHelper.basis1().block()[i],
							modelHelper.basis1().block()[j]);
					SparseElementType tmp = dmrgGeometry_.calcConnectorValue(type,modelHelper.basis1().block()[i],
							dof1,modelHelper.basis1().block()[j],dof2,smax,emin,connectionType);

					if (tmp==0.0) continue;
					
					type = dmrgGeometry_.calcConnectorType(modelHelper.basis1().block()[i],
									modelHelper.basis1().block()[j],modelHelper.basis2().block());	
					if (type==DmrgGeometryType::SystemSystem  ||
						type==DmrgGeometryType::EnvironEnviron) continue; // already included
					flag += (connectionType + 1);
					
					if (lps!=0) {
						lps->isaved.push_back(i);
						lps->jsaved.push_back(j);
						lps->dof1saved.push_back(dof1);
						lps->dof2saved.push_back(dof2);
						lps->typesaved.push_back(type);
						lps->tmpsaved.push_back(tmp);
						lps->connectionsaved.push_back(connectionType);
					} else {
						LinkProductType::calcBond(i,dof1,j,dof2,type,tmp,*matrixBlock,
								modelHelper,connectionType);
					}
					
				}
				return flag;
			}
			
// 			void saveToDisk(VerySparseMatrix<MatrixElementType>& vsm,const std::string& rootname,size_t counter) const
// 			{
// 				typename IoSimple::Out outHandle(rootname+utils::ttos(counter)+".txt",0);
// 				VerySparseMatrix<MatrixElementType> vsmPure(vsm,1e-8);
// 				vsm.clear();
// 				vsmPure.saveToDisk(outHandle);
// 				vsmPure.clear();
// 			}
			
	};     //class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
