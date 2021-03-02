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

	virtual void getMatrix(MatrixType&, SizeType) const = 0;

	virtual SizeType t3OfState(SizeType) const = 0;

	virtual SizeType t8OfState(SizeType) const = 0;

	virtual SizeType size() const = 0;
};

#endif // Su3REPRESENTATION_BASE_H
