/*
Copyright (c) 2009-2012-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

#include "Su2SymmetryGlobals.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "ProgressIndicator.h"
#include "NoPthreads.h"
#include "Sort.h"
#include "Profiling.h"
#include "LinkProductBase.h"
#include "HamiltonianConnection.h"
#include "ParallelHamiltonianConnection.h"

namespace Dmrg {

template<typename ParametersType, typename GeometryType, typename ModelHelperType>
class ModelCommon  {

	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename ModelHelperType::LinkType LinkType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;

public:

	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef LinkProductBase<ModelHelperType, GeometryType> LinkProductBaseType;
	typedef PsimagLite::InputNg<InputCheck>::Readable InputValidatorType;
	typedef typename ModelHelperType::OperatorsType OperatorsType;
	typedef typename ModelHelperType::BlockType Block;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::BasisType MyBasis;
	typedef typename ModelHelperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef HamiltonianConnection<LinkProductBaseType> HamiltonianConnectionType;
	typedef typename HamiltonianConnectionType::CachedHamiltonianLinksType
	CachedHamiltonianLinksType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef typename PsimagLite::Vector<CachedHamiltonianLinksType>::Type
	VectorLinkProductStructType;
	typedef typename HamiltonianConnectionType::VectorSizeType VectorSizeType;
	typedef typename HamiltonianConnectionType::VerySparseMatrixType VerySparseMatrixType;
	typedef ParallelHamiltonianConnection<HamiltonianConnectionType> ParallelHamConnectionType;

	ModelCommon(const ParametersType& params,
	            const GeometryType& geometry,
	            const LinkProductBaseType* lpb)
	    : params_(params),
	      geometry_(geometry),
	      lpb_(lpb),
	      progress_("ModelCommon")
	{
		if (lpb->terms() > geometry.terms()) {
			PsimagLite::String str("ModelCommon: NumberOfTerms must be ");
			err( str + ttos(lpb->terms()) + " in input file for this model\n");
		}

		Su2SymmetryGlobals<RealType>::init(ModelHelperType::isSu2());
		MyBasis::useSu2Symmetry(ModelHelperType::isSu2());
	}

	~ModelCommon()
	{
		delete lpb_;
		lpb_ = 0;
	}

	const LinkProductBaseType& linkProduct() const
	{
		assert(lpb_);
		return *lpb_;
	}

	const ParametersType& params() const { return params_; }

	const GeometryType& geometry() const { return geometry_; }

	/** Let H be the hamiltonian of the  model for basis1 and partition m
	 * consisting of the external product
		 * of basis2 \otimes basis3
		 * This function does x += H*y
		 * The \cppFunction{matrixVectorProduct} function implements the operation $x+=Hy$.
		 * This function
		 * has a default implementation.
		 */
	void matrixVectorProduct(VectorType& x,
	                         const VectorType& y,
	                         const HamiltonianConnectionType& hc) const
	{
		typedef PsimagLite::Parallelizer<ParallelHamConnectionType> ParallelizerType;
		ParallelizerType parallelConnections(PsimagLite::Concurrency::codeSectionParams);

		ParallelHamConnectionType phc(x, y, hc);
		parallelConnections.loopCreate(phc);

		phc.sync();
	}

	/**
		The function \cppFunction{addHamiltonianConnection} implements
		the Hamiltonian connection (e.g. tight-binding links in the case of the Hubbard Model
		or products $S_i\cdot S_j$ in the case of the Heisenberg model) between
		two basis, $basis2$ and $basis3$, in the order of the outer product,
		$basis1={\rm SymmetryOrdering}(basis2\otimes basis3)$. This was
		explained before in Section~\ref{subsec:dmrgBasisWithOperators}.
		This function has a default implementation.
		*/
	void addHamiltonianConnection(SparseMatrixType &matrix,
	                              const LeftRightSuperType& lrs,
	                              RealType currentTime) const
	{
		PsimagLite::Profiling profiling("addHamiltonianConnection", std::cout);
		assert(lrs.super().partition() > 0);
		SizeType total = lrs.super().partition()-1;

		typename PsimagLite::Vector<VerySparseMatrixType*>::Type vvsm(total, 0);
		VectorSizeType nzs(total, 0);

		for (SizeType m = 0; m < total; ++m) {
			SizeType offset = lrs.super().partition(m);
			assert(lrs.super().partition(m + 1) >= offset);
			SizeType bs = lrs.super().partition(m + 1) - offset;

			vvsm[m] = new VerySparseMatrixType(bs, bs);
			VerySparseMatrixType& vsm = *(vvsm[m]);
			HamiltonianConnectionType hc(m, lrs, geometry_, *lpb_, currentTime, 0);

			hc.matrixBond(vsm);
			nzs[m] = vsm.nonZeros();
			if (nzs[m] > 0) continue;
			delete vvsm[m];
			vvsm[m] = 0;
		}

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType permutation(total, 0);
		sort.sort(nzs, permutation);

		typename PsimagLite::Vector<const SparseMatrixType*>::Type vectorOfCrs;

		assert(total == permutation.size());
		for (SizeType i = 0; i < total; ++i) { // loop over new order

			SizeType m = permutation[i]; // get old index from new index

			if (vvsm[m] == 0) continue;

			const VerySparseMatrixType& vsm = *(vvsm[m]);
			SparseMatrixType matrixBlock2;
			matrixBlock2 = vsm;
			delete vvsm[m];
			vvsm[m] = 0;

			SizeType offset = lrs.super().partition(m);
			SparseMatrixType* full = new SparseMatrixType(matrix.rows(),
			                                              matrix.cols(),
			                                              matrixBlock2.nonZeros());
			fromBlockToFull(*full, matrixBlock2, offset);
			vectorOfCrs.push_back(full);
		}

		if (vectorOfCrs.size() == 0) return;

		vectorOfCrs.push_back(&matrix);
		SizeType effectiveTotal = vectorOfCrs.size();

		VectorType ones(effectiveTotal, 1.0);
		SparseMatrixType sumCrs;
		sum(sumCrs, vectorOfCrs, ones);
		vectorOfCrs.pop_back();
		effectiveTotal = vectorOfCrs.size();
		for (SizeType i = 0; i < effectiveTotal; ++i) {
			delete vectorOfCrs[i];
			vectorOfCrs[i] = 0;
		}

		matrix.swap(sumCrs);
	}

	void addConnectionsInNaturalBasis(SparseMatrixType& hmatrix,
	                                  const VectorOperatorType& cm,
	                                  const Block& block,
	                                  bool sysEnvOnly,
	                                  RealType time) const
	{
		if (block.size() != 1)
			err("addConnectionsInNaturalBasis(): unimplemented\n");
	}


	/**
		Returns H, the hamiltonian for basis1 and partition
		$m$ consisting of the external product of basis2$\otimes$basis3
		Note: Used only for debugging purposes
		*/
	void fullHamiltonian(SparseMatrixType& matrix,
	                     const HamiltonianConnectionType& hc) const
	{
		SparseMatrixType matrixBlock;

		//! contribution to Hamiltonian from current system
		hc.modelHelper().calcHamiltonianPart(matrixBlock,true);
		matrix = matrixBlock;

		//! contribution to Hamiltonian from current envirnoment
		hc.modelHelper().calcHamiltonianPart(matrixBlock,false);
		matrix += matrixBlock;

		matrixBlock.clear();

		VerySparseMatrixType vsm(matrix);
		hc.matrixBond(vsm);

		matrix = vsm;
	}

private:

	const ParametersType& params_;
	const GeometryType& geometry_;
	const LinkProductBaseType* lpb_;
	PsimagLite::ProgressIndicator progress_;
}; //class ModelCommon
} // namespace Dmrg
/*@}*/
#endif
