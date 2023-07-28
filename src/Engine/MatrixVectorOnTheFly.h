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
/** \ingroup DMRG */
/*@{*/

/*! \file MatrixVectorOnTheFly.h
 *
 *  A class to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef MATRIX_VECTOR_OTF_H
#define MATRIX_VECTOR_OTF_H

#include "MatrixVectorBase.h"
#include <vector>

namespace Dmrg
{
template <typename ModelType_>
class MatrixVectorOnTheFly : public MatrixVectorBase<ModelType_>
{

	typedef MatrixVectorBase<ModelType_> BaseType;

public:

	typedef ModelType_ ModelType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type value_type;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> FullMatrixType;
	typedef typename ModelType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename ModelHelperType::Aux AuxType;

	MatrixVectorOnTheFly(const ModelType& model,
	    const HamiltonianConnectionType& hc,
	    const AuxType& aux)
	    : model_(model)
	    , hc_(hc)
	    , aux_(aux)
	{
		int maxMatrixRankStored = model.params().maxMatrixRankStored;
		if (hc.modelHelper().size(aux_.m()) > maxMatrixRankStored)
			return;

		model.fullHamiltonian(matrixStored_, hc, aux_);
		assert(isHermitian(matrixStored_, true));
	}

	SizeType rows() const { return hc_.modelHelper().size(aux_.m()); }

	template <typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType& x, SomeVectorType const& y) const
	{
		if (matrixStored_.rows() > 0)
			matrixStored_.matrixVectorProduct(x, y);
		else
			model_.matrixVectorProduct(x, y, hc_, aux_);
	}

	void fullDiag(VectorRealType& eigs, FullMatrixType& fm) const
	{
		int mrs = model_.params().maxMatrixRankStored;
		if (mrs < static_cast<int>(rows())) {
			std::cerr << "Full diag will likely fail, it would need ";
			std::cerr << rows() << " but you gave only " << mrs << "\n";
		}

		BaseType::fullDiag(eigs, fm, matrixStored_, mrs);
	}

private:

	const ModelType& model_;
	const HamiltonianConnectionType& hc_;
	const AuxType& aux_;
	SparseMatrixType matrixStored_;
}; // class MatrixVectorOnTheFly
} // namespace Dmrg

/*@}*/
#endif
