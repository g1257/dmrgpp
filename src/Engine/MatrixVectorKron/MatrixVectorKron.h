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

/*! \file MatrixVectorKron.h
 *
 *  A class to encapsulate the product x+=Hy,
 *  where x and y are vectors and H is the Hamiltonian matrix
 *
 */
#ifndef	MATRIX_VECTOR_KRON_H
#define MATRIX_VECTOR_KRON_H

#include "Vector.h"
#include "InitKron.h"
#include "KronMatrix.h"
#include "MatrixVectorBase.h"

namespace Dmrg {
template<typename ModelType_>
class MatrixVectorKron : public MatrixVectorBase<ModelType_> {

	typedef MatrixVectorBase<ModelType_> BaseType;

	static const bool CHECK_KRON = true;

public:

	typedef ModelType_ ModelType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::RealType RealType;
	typedef typename ModelType::ReflectionSymmetryType ReflectionSymmetryType;
	typedef InitKron<ModelType,ModelHelperType> InitKronType;
	typedef KronMatrix<ModelType,ModelHelperType> KronMatrixType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> FullMatrixType;
	typedef typename SparseMatrixType::value_type value_type;

	MatrixVectorKron(ModelType const *model,
	                 ModelHelperType const *modelHelper,
	                 ReflectionSymmetryType* = 0)
	    : model_(model),initKron_(*model,*modelHelper),kronMatrix_(initKron_)
	{
		int maxMatrixRankStored = model->params().maxMatrixRankStored;
		if (modelHelper->size() > maxMatrixRankStored) return;

		model->fullHamiltonian(matrixStored_,*modelHelper);
		assert(isHermitian(matrixStored_,true));

		checkKron();
	}

	SizeType rank() const { return initKron_.size(); }

	template<typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType &x,SomeVectorType const &y) const
	{
		if (matrixStored_.row() > 0)
			matrixStored_.matrixVectorProduct(x,y);
		else
			kronMatrix_.matrixVectorProduct(x,y);
	}

	void fullDiag(VectorRealType& eigs,FullMatrixType& fm) const
	{
		BaseType::fullDiag(eigs,fm,matrixStored_,model_->params().maxMatrixRankStored);
	}

private:


	void checkKron() const
	{
		if (!CHECK_KRON)
			return;

#ifndef NDEBUG
		return;
#endif

		SizeType n = rank();
		std::cout<<n<<"\n";
		FullMatrixType m(n, n);
		for (SizeType i = 0; i < n; ++i) {
			VectorType e(n, 0.0);
			e[i] = 1.0;
			VectorType ey(n, 0.0);
			kronMatrix_.matrixVectorProduct(ey,e);
			for (SizeType j = 0; j < n; ++j)
				m(i, j) = ey[j];

		}

		std::cout<<"Matrix as thought of by Kron\n";
		std::cout<<m;
		assert(isHermitian(m));
		std::cout<<"Correct matrix\n";
		std::cout<<matrixStored_;
	}

	const ModelType* model_;
	InitKronType initKron_;
	KronMatrixType kronMatrix_;
	SparseMatrixType matrixStored_;
}; // class MatrixVectorKron
} // namespace Dmrg

/*@}*/
#endif

