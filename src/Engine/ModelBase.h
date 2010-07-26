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

/*! \file ModelBase.h
 *
 *  An abstract class to represent the strongly-correlated-electron models that 
 *  can be used with the DmrgSolver
 *
 */
 
#ifndef DMRG_MODEL_BASE
#define DMRG_MODEL_BASE
#include <iostream>

#include "ModelCommon.h"
#include "BasisWithOperators.h"
#include "Su2SymmetryGlobals.h"

namespace Dmrg {
	//! Interface to models for the dmrg solver.
	//! Implement this interface for every new model you wish to write.
	template<typename ModelHelperType,
	typename SparseMatrixType,
 	typename DmrgGeometryType,
  	typename LinkProductType,
	template<typename> class SharedMemoryTemplate>
	class ModelBase  {
		public:
			typedef typename ModelHelperType::OperatorsType OperatorsType;
			typedef typename ModelHelperType::BlockType Block;
			typedef typename ModelHelperType::RealType RealType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef typename ModelHelperType::ReflectionSymmetryType ReflectionSymmetryType;
			typedef typename ModelHelperType::ConcurrencyType ConcurrencyType;
			typedef Basis<RealType,SparseMatrixType> MyBasis;
			typedef BasisWithOperators<OperatorsType,ConcurrencyType> MyBasisWithOperators;
			typedef ModelCommon<ModelHelperType,SparseMatrixType,DmrgGeometryType,
   					LinkProductType,SharedMemoryTemplate> ModelCommonType;
			typedef DmrgGeometryType GeometryType;
			typedef typename ModelCommonType::SharedMemoryType SharedMemoryType;
			
			ModelBase(int DEGREES_OF_FREEDOM,const DmrgGeometryType& dmrgGeometry) :
					modelCommon_(DEGREES_OF_FREEDOM,dmrgGeometry)
			{
				Su2SymmetryGlobals<RealType>::init(ModelHelperType::isSu2());
				MyBasis::useSu2Symmetry(ModelHelperType::isSu2());
			}  

			//! return internal degrees of freedom for this model
			int dof() const;

			//! return the maximum system size
			int absoluteSystemSize() const;

			//! Let H be the hamiltonian of the  model for basis1 and partition m consisting of the external product
			//! of basis2 \otimes basis3
			//! This function does x += H*y
			template<typename SomeVectorType>
			void matrixVectorProduct(SomeVectorType &x,SomeVectorType const &y,ModelHelperType const &modelHelper) const
			{
				modelCommon_.matrixVectorProduct(x,y,modelHelper);
			}

			//! site of the hamiltonian matrix
			int getSize(ModelHelperType const &modelHelper) const; 

			//! electronic density
			double density() const;

			void addHamiltonianConnection(SparseMatrixType &matrix,MyBasis const &basis1,MyBasisWithOperators const &basis2,
				MyBasisWithOperators const &basis3,size_t dof,size_t nOrbitals) const
			{	
				int bs,offset;
				SparseMatrixType matrixBlock;

				for (size_t m=0;m<basis1.partition()-1;m++) {
					offset =basis1.partition(m);
					bs = basis1.partition(m+1)-offset;
					matrixBlock.makeDiagonal(bs);
					ModelHelperType modelHelper(m,basis1,basis2,basis3,nOrbitals);
					modelCommon_.addHamiltonianConnection(matrixBlock,modelHelper);
					sumBlock(matrix,matrixBlock,offset);
				}
			}

			//! Add Hamiltonian connection between basis2 and basis3 in the orderof basis1 for symmetry block m 
			void addHamiltonianConnection(VerySparseMatrix<SparseElementType>& matrix,const ModelHelperType& modelHelper) const
			{
				modelCommon_.addHamiltonianConnection(matrix,modelHelper);
			}

			//! Let H_m be the Hamiltonian connection between basis2 and basis3 in the orderof basis1 for block m 
			//! Then this function does x+= H_m *y
			void hamiltonianConnectionProduct(std::vector<SparseElementType> &x,std::vector<SparseElementType> const &y,
				ModelHelperType const &modelHelper) const
			{
				return modelCommon_.hamiltonianConnectionProduct(x,y,modelHelper);
			}

			//! find  operator matrices for in the natural basis and quantum numbers
			void setNaturalBasis(std::vector<SparseMatrixType> &operatorMatrices,SparseMatrixType &hamiltonian,
					     std::vector<int> &q,std::vector<size_t> &electrons,Block const &block) const;
			
			//! set operator matrices for sites in block
			void setOperatorMatrices(std::vector<SparseMatrixType> &operatorMatrix,Block const &block) const;

			//! split into S X Y and E
			void setBlocksOfSites(Block &S,std::vector<Block> &X,std::vector<Block> &Y,Block &E) const; 

			//! print model or model parameters
			void print(std::ostream& os) const;

			//! Given a vector of quantum numbers, return the effective quantum number 
			//! (see paper for definition of "effective quantum number")
			//! THis is used to run a single symmetry sector as specified in the input file
			int getTargetSector(const std::vector<int>& quantumNumbers) const;
			
			//! Return H, the hamiltonian of the FeAs model for basis1 and partition m consisting of the external product
			//! of basis2 \otimes basis3
			//! Note: Used only for debugging purposes
			void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
			{
				modelCommon_.fullHamiltonian(matrix,modelHelper);
			}
	
	
		private:
			ModelCommonType modelCommon_;
	};     //class ModelBase
} // namespace Dmrg
/*@}*/
#endif
