#ifndef SU3REPRESENTATIONP1_H
#define SU3REPRESENTATIONP1_H
#include "PsimagLite.h"
#include "Matrix.h"
#include "Su3RepresentationBase.h"

template<typename ComplexOrRealType>
class Su3RepresentationP1 : public Su3RepresentationBase<ComplexOrRealType> {

public:

	//typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	void getMatrix(MatrixType& m, SizeType n) const
	{
		throw PsimagLite::RuntimeError("getMatrix: not implemented yet\n");
	}

	SizeType t3OfState(SizeType) const
	{
		throw PsimagLite::RuntimeError("t3OfState: not implemented yet\n");
	}

	SizeType t8OfState(SizeType) const
	{
		throw PsimagLite::RuntimeError("t8OfState: not implemented yet\n");
	}

	SizeType size() const { return 3; }
};

#endif // SU3REPRESENTATIONP1_H
