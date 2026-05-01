/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef DEFAULT_SYMM_H
#define DEFAULT_SYMM_H
#include "CrsMatrix.h"
#include "LanczosGlobals.h"
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "Vector.h"
#include <iostream>

namespace LanczosPlusPlus {

template <typename GeometryType_, typename BasisType> class DefaultSymmetry {

	typedef typename GeometryType_::ComplexOrRealType          ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef LanczosGlobals::WordType                           WordType;
	typedef PsimagLite::Matrix<ComplexOrRealType>              MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type        VectorRealType;

public:

	typedef GeometryType_                                        GeometryType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType>             SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<VectorType>::Type        VectorVectorType;

	DefaultSymmetry(const BasisType&, const GeometryType&, PsimagLite::String options)
	    : printMatrix_(options.find("printmatrix") != PsimagLite::String::npos)
	    , dumpMatrix_(options.find("dumpmatrix") != PsimagLite::String::npos)
	{ }

	template <typename SomeModelType>
	void init(const SomeModelType& model, const BasisType& basis)
	{
		model.setupHamiltonian(matrixStored_, basis);
		int nrows = matrixStored_.rows();
		if (printMatrix_ && nrows > 512)
			throw PsimagLite::RuntimeError("printMatrix: too big\n");

		if (printMatrix_ || dumpMatrix_) {
			std::cout << "#LanczosPlusPlus: Basis for matrix\n";

			if (printMatrix_) {
				basis.print(std::cout, BasisType::PRINT_BINARY);
				std::cout << "DenseMatrix\n";
				std::cout << matrixStored_.toDense();
			}

			basis.print(std::cout, BasisType::PRINT_DECIMAL);
			model.printOperators(std::cout);

			assert(isHermitian(matrixStored_, true));
			PsimagLite::Matrix<ComplexOrRealType> matrixCopy;
			VectorRealType                        eigs;
			fullDiag(eigs, matrixCopy);
		}
	}

	void fullDiag(VectorRealType& eigs, MatrixType& fm) const
	{
		if (matrixStored_.rows() > 4900)
			throw PsimagLite::RuntimeError("fullDiag too big\n");

		fm = matrixStored_.toDense();
		diag(fm, eigs, 'V');

		if (printMatrix_ || dumpMatrix_) {
			std::cout << "#Eigenvalues\n";
			std::cout << eigs;
			std::cout << "#Eigenvectors\n";
			std::cout << fm;
		}
	}

	void transformMatrix(typename PsimagLite::Vector<SparseMatrixType>::Type&,
	                     const SparseMatrixType&) const
	{
		throw std::runtime_error("DefaultSymmetry: cannot call transformMatrix\n");
	}

	void transform(VectorVectorType&, SizeType) { }

	SizeType sectors() const { return 1; }

	void setPointer(SizeType) { }

	PsimagLite::String name() const { return "default"; }

	SizeType rows() const { return matrixStored_.rows(); }

	template <typename SomeVectorType>
	void matrixVectorProduct(SomeVectorType& x, SomeVectorType const& y) const
	{
		return matrixStored_.matrixVectorProduct(x, y);
	}

private:

	SparseMatrixType matrixStored_;
	bool             printMatrix_;
	bool             dumpMatrix_;
}; // class DefaultSymmetry
} // namespace Dmrg

#endif // DEFAULT_SYMM_H
