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

/*! \file ModelCommon.h
 *
 *  Common functions for most models
 *
 */
 
#ifndef MODEL_COMMON_H
#define MODEL_COMMON_H

#include "VerySparseMatrix.h"
#include "IoSimple.h"
#include "HamiltonianConnection.h"

namespace Dmrg {
	//! Common functions for various models
	
	template<typename ModelHelperType,
	typename SparseMatrixType, // <-- kill this template, deduced from ModelHelper FIXME
 	typename DmrgGeometryType,
	typename LinkProductType,
	template<typename> class SharedMemoryTemplate>
	class ModelCommon  {
		public:	
			typedef typename ModelHelperType::RealType RealType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef VerySparseMatrix<SparseElementType> VerySparseMatrixType;
			typedef HamiltonianConnection<DmrgGeometryType,ModelHelperType,
					LinkProductType> HamiltonianConnectionType;
			typedef typename HamiltonianConnectionType::LinkProductStructType LinkProductStructType;
			typedef SharedMemoryTemplate<HamiltonianConnectionType> SharedMemoryType;
			
			ModelCommon(const DmrgGeometryType& dmrgGeometry)
			: dmrgGeometry_(dmrgGeometry)
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
			template<typename SomeModelHelperType>
			void addHamiltonianConnection(
					SparseMatrixType& matrix,
					const SomeModelHelperType& modelHelper) const
			{
				size_t n = matrix.rank();
				VerySparseMatrixType vsm(n);
				addHamiltonianConnection(vsm,modelHelper);
				SparseMatrixType matrixBlock;
				matrixBlock = vsm;
				matrix += matrixBlock;
			}
			
			size_t maxConnections() const
			{
				return dmrgGeometry_.maxConnections();
			}

			//! Let H_m be the Hamiltonian connection between basis2 and basis3 in the orderof basis1 for block m 
			//! Then this function does x+= H_m *y
			void hamiltonianConnectionProduct(
					std::vector<SparseElementType> &x,
					const std::vector<SparseElementType>& y,
					const ModelHelperType& modelHelper) const
			{
				size_t n=modelHelper.leftRightSuper().super().block().size();

				//SparseMatrixType matrix;
				size_t maxSize = maxConnections() * 4 * 16;
				maxSize *= maxSize;
				
				static LinkProductStructType lps(maxSize);
				HamiltonianConnectionType hc(dmrgGeometry_,modelHelper,&lps,&x,&y);
				
				size_t total = 0;
				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
						hc.compute(i,j,0,&lps,total);
					}
				}

				SharedMemoryType pthreads;
				pthreads.loopCreate(total,hc);

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
// 				PsimagLite::Matrix<SparseElementType> fm;
// 				crsMatrixToFullMatrix(fm,matrix);
// 				std::vector<RealType> e(matrix.rank());
// 				printNonZero(fm,std::cerr);
// 				PsimagLite::diag(fm,e,'N');
// 				std::cerr<<"eSystem[0]="<<e[0]<<"\n";
	
				//! contribution to Hamiltonian from current envirnoment
				modelHelper.calcHamiltonianPart(matrixBlock,false);
				matrix += matrixBlock;
// 				crsMatrixToFullMatrix(fm,matrixBlock);
// 				printNonZero(fm,std::cerr);
// 				PsimagLite::diag(fm,e,'N');
// 				std::cerr<<"eEnv[0]="<<e[0]<<"\n";
				matrixBlock.clear();
	
				VerySparseMatrixType vsm(matrix);
				addHamiltonianConnection(vsm,modelHelper);
	
				matrix = vsm;
			}

		private:

			//! Add Hamiltonian connection between basis2 and basis3 in
			// the orderof basis1 for symmetry block m
			template<typename SomeModelHelperType>
			void addHamiltonianConnection(VerySparseMatrixType& matrix,
					const SomeModelHelperType& modelHelper) const
			{
				size_t n=modelHelper.leftRightSuper().sites();
				size_t matrixRank = matrix.rank();
				VerySparseMatrixType matrix2(matrixRank);
				typedef HamiltonianConnection<
						DmrgGeometryType,
						SomeModelHelperType,
						LinkProductType> SomeHamiltonianConnectionType;
				SomeHamiltonianConnectionType hc(dmrgGeometry_,modelHelper);

				size_t total = 0;
				for (size_t i=0;i<n;i++) {
					for (size_t j=0;j<n;j++) {
						SparseMatrixType matrixBlock(matrixRank,matrixRank);
						if (!hc.compute(i,j,&matrixBlock,0,total)) continue;
						VerySparseMatrixType vsm(matrixBlock);
						matrix2+=vsm;
					}
				}
				matrix += matrix2;
			}

			const DmrgGeometryType& dmrgGeometry_;
			
	};     //class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
