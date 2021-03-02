#ifndef Su3REPRESENTATION_BASE_H
#define Su3REPRESENTATION_BASE_H
#include "PsimagLite.h"
#include "Matrix.h"

template<typename ComplexOrRealType>
class Su3RepresentationBase {

public:

	//typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	virtual ~Su3RepresentationBase() {}

	virtual void getMatrix(MatrixType& m, SizeType n) const = 0;

	virtual SizeType size() const = 0;
};

#endif // Su3REPRESENTATION_BASE_H
