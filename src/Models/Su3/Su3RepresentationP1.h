#ifndef SU3REPRESENTATIONP1_H
#define SU3REPRESENTATIONP1_H
#include "PsimagLite.h"
#include "Matrix.h"
#include "Su3RepresentationBase.h"

template<typename ComplexOrRealType, bool>
class Su3RepresentationP1 : public Su3RepresentationBase<ComplexOrRealType> {

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	void getMatrix(MatrixType& m, SizeType n) const
	{
		m.resize(3, 3);

		// Tplus1
		if (n == 0) {
			m(0, 1) = 0.5;
			return;
		}

		// Tplus2
		if (n == 1) {
			m(0, 2) = 0.5;
			return;
		}

		// Tplus3
		if (n == 2) {
			m(1, 2) = 0.5;
			return;
		}

		// T3
		if (n == 3) {
			m(0, 0) = 0.5;
			m(1, 1) = -0.5;
			return;
		}

		// T8
		assert(n == 4);
		static const ComplexOrRealType oneOverSqrt3 = 1.0/sqrt(3.0);
		m(0, 0) = m(1, 1) = 0.5*oneOverSqrt3;
		m(2, 2) = -oneOverSqrt3;
	}

	// diagonal 1 -1 0  ==> 2 0 1
	SizeType t3OfState(SizeType ind) const
	{
		if (ind == 0) return 2;
		if (ind == 1) return 0;
		assert(ind == 2);
		return 1;
	}

	// diagonal 1 1 -2 ==> 3 3 0
	SizeType t8OfState(SizeType ind) const
	{
		assert(ind < 3);
		return (ind < 2) ? 3 : 0;
	}

	SizeType size() const { return 3; }
};

template<typename ComplexOrRealType>
class Su3RepresentationP1<ComplexOrRealType, true>
        : public Su3RepresentationBase<ComplexOrRealType> {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

	void getMatrix(MatrixType& m, SizeType n) const
	{
		m.resize(3, 3);

		if (n == 0) {
			m(0, 1) = m(1, 0) = 0.5;
			return;
		}

		if (n == 1) {
			m(0, 1) = std::complex<RealType>(0, -0.5);
			m(1, 0) = -m(0, 1);
			return;
		}

		if (n == 2) {
			m(0, 0) = 0.5;
			m(1, 1) = -0.5;
			return;
		}

		if (n == 3) {
			m(0, 2) = m(2, 0) = 0.5;
			return;
		}

		if (n == 4) {
			m(0, 2) = std::complex<RealType>(0, -0.5);
			m(2, 0) = -m(0, 2);
			return;
		}

		if (n == 5) {
			m(1, 2) = m(2, 1) = 0.5;
			return;
		}

		if (n == 6) {
			m(1, 2) = std::complex<RealType>(0, -0.5);
			m(2, 1) = -m(1, 2);
			return;
		}

		assert(n == 7);
		static const RealType oneOverSqrt3 = 1.0/sqrt(3.0);
		m(0, 0) = m(1, 1) = 0.5*oneOverSqrt3;
		m(2, 2) = -oneOverSqrt3;
	}

	// diagonal 1 -1 0  ==> 2 0 1
	SizeType t3OfState(SizeType ind) const
	{
		if (ind == 0) return 2;
		if (ind == 1) return 0;
		assert(ind == 2);
		return 1;
	}

	// diagonal 1 1 -2 ==> 3 3 0
	SizeType t8OfState(SizeType ind) const
	{
		assert(ind < 3);
		return (ind < 2) ? 3 : 0;
	}

	SizeType size() const { return 3; }
};

#endif // SU3REPRESENTATIONP1_H
