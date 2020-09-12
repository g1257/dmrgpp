#ifndef EXACTGREENFUNCTION_H
#define EXACTGREENFUNCTION_H
#include "Vector.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include "CrsMatrix.h"

namespace Dmft {

template<typename ComplexOrRealType>
class ExactGreenFunction {

public:

	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorComplexType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;

	const VectorComplexType& operator()() const { return gimp_; }

	// <gs|c'(iwn-Hbar)^{-1}c|gs> + <gs|c(iwn+Hbar)^{-1}c'|gs>
	void operator()(RealType energy, const VectorComplexType& gs, const SparseMatrixType& sparse)
	{
		PsimagLite::Matrix<ComplexOrRealType> mdense = sparse.toDense();
		err("ExactGreenFunction not implemented yet");
	}

private:

	VectorComplexType gimp_;

};
}
#endif // EXACTGREENFUNCTION_H
