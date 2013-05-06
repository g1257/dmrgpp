/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file ModelBase.h
 *
 *  An abstract class to represent the strongly-correlated-electron models that 
 *  can be used with the DmrgSolver
 *
 */
 
#ifndef DMRG_MODEL_BASE
#define DMRG_MODEL_BASE
#include <iostream>

#include "VerySparseMatrix.h"
#include "IoSimple.h"
#include "HamiltonianConnection.h"
#include "Su2SymmetryGlobals.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "ProgressIndicator.h"

namespace Dmrg {
	
	/**
	A sample SCE model, the one-band Hubbard model,
	\[
	\sum_{i,j,\sigma}t_{ij}c^\dagger_{i\sigma} c_{j\sigma}
	+\sum_i U_i n_{i\uparrow}n_{i\downarrow} + \sum_{i,\sigma}V_i n_{i\sigma},
	\]
	is implemented in class \cppClass{ModelHubbardOneBand}. 
	A sample \cppClass{HeisenbergSpinOneHalf} is also included for the Heisenberg model 
	$\sum_{ij}J_{ij}\vec{S}_i\cdot\vec{S}_j$.
	These models  inherit from the abstract class \cppClass{!PTEX_THISCLASS}. 
	To implement other SCE models one has to implement the functions prototyped 
	in this abstract class. Note that there are default implementations for
	some of these functions; they delegate to the \cppClass{ModelCommon} class.
	Interface (functions in \cppClass{!PTEX_THISCLASS}) are the following.
	*/
	template<typename ModelHelperType,
	typename SparseMatrixType,
 	typename DmrgGeometryType,
  	typename LinkProductType,
	template<typename> class ParallelConnectionsTemplate>
	class ModelBase  {

		typedef typename SparseMatrixType::value_type SparseElementType;
		typedef VerySparseMatrix<SparseElementType> VerySparseMatrixType;

	public:

		typedef PsimagLite::InputNg<InputCheck>::Readable InputValidatorType;
		typedef typename ModelHelperType::OperatorsType OperatorsType;
		typedef typename ModelHelperType::BlockType Block;
		typedef typename ModelHelperType::RealType RealType;
		typedef typename ModelHelperType::ConcurrencyType
				ConcurrencyType;
		typedef typename ModelHelperType::BasisType MyBasis;
		typedef typename ModelHelperType::BasisWithOperatorsType
				BasisWithOperatorsType;
		typedef DmrgGeometryType GeometryType;
		typedef HamiltonianConnection<DmrgGeometryType,ModelHelperType,LinkProductType> HamiltonianConnectionType;
		typedef ParallelConnectionsTemplate<HamiltonianConnectionType> ParallelConnectionsType;
		typedef typename HamiltonianConnectionType::LinkProductStructType LinkProductStructType;
		typedef typename ModelHelperType::LeftRightSuperType
				LeftRightSuperType;
		typedef typename OperatorsType::OperatorType OperatorType;
		typedef typename MyBasis::BasisDataType BasisDataType;

		ModelBase(const DmrgGeometryType& geometry,ConcurrencyType& concurrency)
		: dmrgGeometry_(geometry),concurrency_(concurrency),progress_("ModelBase",concurrency.rank())
		{
			Su2SymmetryGlobals<RealType>::init(ModelHelperType::isSu2());
			MyBasis::useSu2Symmetry(ModelHelperType::isSu2());
		}

		/** Let H be the hamiltonian of the  model for basis1 and partition m consisting of the external product
		 * of basis2 \otimes basis3
		 * This function does x += H*y
		 * The \\cppFunction{matrixVectorProduct} function implements the operation $x+=Hy$. This function
		 * has a default implementation.
		 */
		void matrixVectorProduct(typename PsimagLite::Vector<RealType>::Type& x,
		                         const typename PsimagLite::Vector<RealType>::Type& y,
		                         ModelHelperType const &modelHelper) const
		{
			//! contribution to Hamiltonian from current system
			modelHelper.hamiltonianLeftProduct(x,y);
			//! contribution to Hamiltonian from current envirnoment
			modelHelper.hamiltonianRightProduct(x,y);
			//! contribution to Hamiltonian from connection system-environment
			hamiltonianConnectionProduct(x,y,modelHelper);
		}

		void matrixVectorProduct(typename PsimagLite::Vector<std::complex<RealType> >::Type& x,
		                         const typename PsimagLite::Vector<std::complex<RealType> >::Type& y,
		                         ModelHelperType const &modelHelper) const
		{
			//! contribution to Hamiltonian from current system
			modelHelper.hamiltonianLeftProduct(x,y);
			//! contribution to Hamiltonian from current envirnoment
			modelHelper.hamiltonianRightProduct(x,y);
			//! contribution to Hamiltonian from connection system-environment
			hamiltonianConnectionProduct(x,y,modelHelper);
		}

		/**
		The function \cppFunction{addHamiltonianConnection} implements
		the Hamiltonian connection (e.g. tight-binding links in the case of the Hubbard Model
		or products $S_i\cdot S_j$ in the case of the Heisenberg model) between 
		two basis, $basis2$ and $basis3$, in the order of the outer product,
		$basis1={\\rm SymmetryOrdering}(basis2\otimes basis3)$. This was
		explained before in Section~\\ref{subsec:dmrgBasisWithOperators}.
		This function has a default implementation.
		*/
		void addHamiltonianConnection(SparseMatrixType &matrix,const LeftRightSuperType& lrs) const
		{
			int bs,offset;
			SparseMatrixType matrixBlock;

			for (size_t m=0;m<lrs.super().partition()-1;m++) {
				offset =lrs.super().partition(m);
				bs = lrs.super().partition(m+1)-offset;
				matrixBlock.makeDiagonal(bs);
				ModelHelperType modelHelper(m,lrs);

				VerySparseMatrixType vsm(matrixBlock.row());
				addHamiltonianConnection(vsm,modelHelper);
				SparseMatrixType matrixBlock2;
				matrixBlock2 = vsm;
				matrixBlock += matrixBlock2;

				sumBlock(matrix,matrixBlock,offset);
			}
		}

		size_t maxConnections() const
		{
			return dmrgGeometry_.maxConnections();
		}
		
		/**
		Let $H_m$ be the Hamiltonian connection between basis2 and basis3 in 
		the orderof basis1 for block $m$. Then this function does $x+= H_m *y$
		*/
		void hamiltonianConnectionProduct(typename PsimagLite::Vector<SparseElementType>::Type& x,
		                                  const typename PsimagLite::Vector<SparseElementType>::Type& y,
		                                  ModelHelperType const &modelHelper) const
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

			ParallelConnectionsType parallelConnections;
			parallelConnections.loopCreate(total,hc,concurrency_);
			hc.sync(parallelConnections,concurrency_);
		}

		/**
		Returns H, the hamiltonian for basis1 and partition 
		$m$ consisting of the external product of basis2$\\otimes$basis3
		Note: Used only for debugging purposes
		*/
		void fullHamiltonian(SparseMatrixType& matrix,const ModelHelperType& modelHelper) const
		{
			SparseMatrixType matrixBlock;

			//! contribution to Hamiltonian from current system
			modelHelper.calcHamiltonianPart(matrixBlock,true);
			matrix = matrixBlock;

			//! contribution to Hamiltonian from current envirnoment
			modelHelper.calcHamiltonianPart(matrixBlock,false);
			matrix += matrixBlock;

			matrixBlock.clear();

			VerySparseMatrixType vsm(matrix);
			addHamiltonianConnection(vsm,modelHelper);

			matrix = vsm;
		}

	private:

		//! Add Hamiltonian connection between basis2 and basis3 in the orderof basis1 for symmetry block m
		void addHamiltonianConnection(VerySparseMatrix<SparseElementType>& matrix,
					      const ModelHelperType& modelHelper) const
		{
			size_t n=modelHelper.leftRightSuper().sites();
			size_t matrixRank = matrix.rank();
			VerySparseMatrixType matrix2(matrixRank);
			typedef HamiltonianConnection<
					DmrgGeometryType,
					ModelHelperType,
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
		ConcurrencyType& concurrency_;
		PsimagLite::ProgressIndicator progress_;
	};     //class ModelBase
} // namespace Dmrg
/*@}*/
#endif
