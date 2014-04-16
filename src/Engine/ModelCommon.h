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

/*! \file ModelCommon.h
 *
 *  An abstract class to represent the strongly-correlated-electron models that
 *  can be used with the DmrgSolver
 *
 */

#ifndef MODEL_COMMON_H
#define MODEL_COMMON_H
#include <iostream>

#include "VerySparseMatrix.h"
#include "IoSimple.h"
#include "HamiltonianConnection.h"
#include "Su2SymmetryGlobals.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "ProgressIndicator.h"

namespace Dmrg {


template<typename ModelBaseType,typename LinkProductType>
class ModelCommon : public ModelBaseType::ModelCommonBaseType {

	typedef typename ModelBaseType::ModelCommonBaseType ModelCommonBaseType;
	typedef typename ModelBaseType::ModelHelperType ModelHelperType;
	typedef typename ModelBaseType::GeometryType GeometryType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef VerySparseMatrix<SparseElementType> VerySparseMatrixType;
	typedef typename ModelHelperType::LinkType LinkType;

public:

	typedef PsimagLite::InputNg<InputCheck>::Readable InputValidatorType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::BlockType Block;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::BasisType MyBasis;
	typedef typename ModelHelperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef HamiltonianConnection<GeometryType,ModelHelperType,LinkProductType> HamiltonianConnectionType;
	typedef typename HamiltonianConnectionType::LinkProductStructType LinkProductStructType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename MyBasis::BasisDataType BasisDataType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename ModelBaseType::SolverParamsType SolverParamsType;
	typedef typename PsimagLite::Vector<LinkProductStructType>::Type VectorLinkProductStructType;

	ModelCommon(const SolverParamsType& params,const GeometryType& geometry)
	    : ModelCommonBaseType(params,geometry),progress_("ModelCommon")
	{
		if (LinkProductType::terms() > this->geometry().terms()) {
			PsimagLite::String str("ModelCommon: NumberOfTerms must be ");
			str += ttos(LinkProductType::terms()) + " in input file for this model\n";
			throw PsimagLite::RuntimeError(str);
		}

		Su2SymmetryGlobals<RealType>::init(ModelHelperType::isSu2());
		MyBasis::useSu2Symmetry(ModelHelperType::isSu2());

		SizeType total = 1;
		PsimagLite::String options = this->params().options;
		bool cTridiag = (options.find("concurrenttridiag") != PsimagLite::String::npos);
		if (cTridiag) total = PsimagLite::Concurrency::npthreads;

		SizeType maxSize = this->geometry().maxConnections() * 4 * 16;
		maxSize *= maxSize;
		vlps_.resize(total,maxSize);

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
		//! contribution to Hamiltonian from connection system-environment
		hamiltonianConnectionProduct(x,y,modelHelper);
		//! contribution to Hamiltonian from current system
		modelHelper.hamiltonianLeftProduct(x,y);
		//! contribution to Hamiltonian from current envirnoment
		modelHelper.hamiltonianRightProduct(x,y);
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

		for (SizeType m=0;m<lrs.super().partition()-1;m++) {
			offset =lrs.super().partition(m);
			bs = lrs.super().partition(m+1)-offset;
			matrixBlock.makeDiagonal(bs);
			SizeType threadId = 0;
			ModelHelperType modelHelper(m,lrs,threadId);

			VerySparseMatrixType vsm(matrixBlock.row());
			addHamiltonianConnection(vsm,modelHelper);
			SparseMatrixType matrixBlock2;
			matrixBlock2 = vsm;
			matrixBlock += matrixBlock2;

			sumBlock(matrix,matrixBlock,offset);
		}
	}

	SizeType maxConnections() const
	{
		return this->geometry().maxConnections();
	}

	/**
		Let $H_m$ be the Hamiltonian connection between basis2 and basis3 in
		the orderof basis1 for block $m$. Then this function does $x+= H_m *y$
		*/
	void hamiltonianConnectionProduct(typename PsimagLite::Vector<SparseElementType>::Type& x,
	                                  const typename PsimagLite::Vector<SparseElementType>::Type& y,
	                                  ModelHelperType const &modelHelper) const
	{
		SizeType threadId = modelHelper.threadId();

		if (threadId >= vlps_.size())
			throw PsimagLite::RuntimeError("hamiltonianConnectionProduct\n");

		LinkProductStructType& lps = vlps_[threadId];
		HamiltonianConnectionType hc(this->geometry(),modelHelper,&lps,&x,&y);

		SizeType n=modelHelper.leftRightSuper().super().block().size();
		SizeType total = 0;
		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<n;j++) {
				hc.compute(i,j,0,&lps,total);
			}
		}

		PsimagLite::String options = this->params().options;
		bool cTridiag = (options.find("concurrenttridiag") != PsimagLite::String::npos);

		if (cTridiag) {
			typedef PsimagLite::NoPthreads<HamiltonianConnectionType> ParallelizerType;
			ParallelizerType parallelConnections(1,0);
			parallelConnections.loopCreate(total,hc);
		} else {
			typedef PsimagLite::Parallelizer<HamiltonianConnectionType> ParallelizerType;
			ParallelizerType parallelConnections(PsimagLite::Concurrency::npthreads,
			                                     PsimagLite::MPI::COMM_WORLD);
			parallelConnections.loopCreate(total,hc);
		}

		hc.sync();

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

	void addConnectionsInNaturalBasis(SparseMatrixType& hmatrix,
	                                  const VectorOperatorType& cm,
	                                  const Block& block,
	                                  bool sysEnvOnly) const
	{
		SizeType n = block.size();
		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<n;j++) {
				addConnectionsInNaturalBasis(hmatrix,i,j,cm,block,sysEnvOnly);
			}
		}
	}

	SizeType getLinkProductStruct(LinkProductStructType** lps,
	                              const ModelHelperType& modelHelper) const
	{
		SizeType n=modelHelper.leftRightSuper().super().block().size();

		SizeType maxSize = vlps_[0].isaved.size();
		*lps = new LinkProductStructType(maxSize);

		typename PsimagLite::Vector<SparseElementType>::Type x,y; // bogus

		HamiltonianConnectionType hc(this->geometry(),modelHelper,*lps,&x,&y);
		SizeType total = 0;
		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<n;j++) {
				hc.compute(i,j,0,*lps,total);
			}
		}
		return total;
	}

	LinkType getConnection(const SparseMatrixType** A,
	                       const SparseMatrixType** B,
	                       SizeType ix,
	                       const LinkProductStructType& lps,
	                       const ModelHelperType& modelHelper) const
	{
		typename PsimagLite::Vector<SparseElementType>::Type x,y; // bogus
		HamiltonianConnectionType hc(this->geometry(),modelHelper,&lps,&x,&y);
		SizeType i =0, j = 0, type = 0,term = 0, dofs =0;
		SparseElementType tmp = 0.0;
		typename HamiltonianConnectionType::AdditionalDataType additionalData;
		hc.prepare(ix,i,j,type,tmp,term,dofs,additionalData);
		LinkType link2 = hc.getKron(A,B,i,j,type,tmp,term,dofs,additionalData);
		return link2;
	}

private:

	void addConnectionsInNaturalBasis(SparseMatrixType& hmatrix,
	                                  SizeType i,
	                                  SizeType j,
	                                  const VectorOperatorType& cm,
	                                  const Block& block,
	                                  bool sysEnvOnly) const
	{
		SizeType middle = static_cast<SizeType>(block.size()/2);
		if (middle == 0 && sysEnvOnly) return;

		assert(middle < block.size());
		middle = block[middle];

		SizeType ind = block[i];
		SizeType jnd = block[j];

		if (!this->geometry().connected(0,1,ind,jnd)) return;
		if (sysEnvOnly && ind < middle && jnd < middle) return;
		if (sysEnvOnly && ind >= middle && jnd >= middle) return;

		SizeType type = 0;
		SizeType offset = cm.size()/block.size();

		typename GeometryType::AdditionalDataType additionalData;

		for (SizeType term=0;term<this->geometry().terms();term++) {
			this->geometry().fillAdditionalData(additionalData,term,ind,jnd);
			SizeType dofsTotal = LinkProductType::dofs(term,additionalData);
			for (SizeType dofs=0;dofs<dofsTotal;dofs++) {
				std::pair<SizeType,SizeType> edofs = LinkProductType::connectorDofs(term,dofs,additionalData);
				SparseElementType tmp = this->geometry()(ind,edofs.first,jnd,edofs.second,term);

				if (tmp==static_cast<RealType>(0.0)) continue;

				std::pair<SizeType,SizeType> ops;
				std::pair<char,char> mods('N','C');
				SizeType fermionOrBoson=ProgramGlobals::FERMION,angularMomentum=0,category=0;
				RealType angularFactor=0;
				bool isSu2 = ModelHelperType::isSu2();
				SparseElementType value = tmp;
				LinkProductType::valueModifier(value,term,dofs,isSu2,additionalData);
				LinkProductType::setLinkData(term,dofs,isSu2,fermionOrBoson,ops,mods,angularMomentum,angularFactor,category,additionalData);
				typename ModelHelperType::LinkType link2(i,j,type, value,dofs,fermionOrBoson,ops,mods,angularMomentum,angularFactor,category);

				const SparseMatrixType& A = cm[link2.ops.first+i*offset].data;

				const SparseMatrixType& B = cm[link2.ops.second+j*offset].data;

				if (ind > jnd) tmp = std::conj(tmp);

				hmatrix += tmp * (transposeOrNot(B,link2.mods.second) *
				                  transposeOrNot(A,link2.mods.first));
			}
		}
	}

	//! Add Hamiltonian connection between basis2 and basis3 in the orderof basis1 for symmetry block m
	void addHamiltonianConnection(VerySparseMatrix<SparseElementType>& matrix,
	                              const ModelHelperType& modelHelper) const
	{
		SizeType n=modelHelper.leftRightSuper().sites();
		SizeType matrixRank = matrix.rank();
		VerySparseMatrixType matrix2(matrixRank);
		typedef HamiltonianConnection<
		        GeometryType,
		        ModelHelperType,
		        LinkProductType> SomeHamiltonianConnectionType;
		SomeHamiltonianConnectionType hc(this->geometry(),modelHelper);

		SizeType total = 0;
		for (SizeType i=0;i<n;i++) {
			for (SizeType j=0;j<n;j++) {
				SparseMatrixType matrixBlock(matrixRank,matrixRank);
				if (!hc.compute(i,j,&matrixBlock,0,total)) continue;
				VerySparseMatrixType vsm(matrixBlock);
				matrix2+=vsm;
			}
		}
		matrix += matrix2;
	}

	SparseMatrixType transposeOrNot(const SparseMatrixType& A,char mod) const
	{
		if (mod == 'C' || mod == 'c') {
			SparseMatrixType Ac;
			transposeConjugate(Ac,A);
			return Ac;
		}
		return A;
	}

	PsimagLite::ProgressIndicator progress_;
	mutable VectorLinkProductStructType vlps_;
};     //class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
