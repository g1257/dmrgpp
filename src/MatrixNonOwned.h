#ifndef DMRG_MATRIX_NON_OWNED_H
#define DMRG_MATRIX_NON_OWNED_H
#include "Matrix.h"

namespace PsimagLite {

template<typename T>
struct VectorConstOrNot {
	typedef typename Vector<T>::Type Type;
};

template<typename T>
struct VectorConstOrNot<const T> {
	typedef const typename Vector<T>::Type Type;
};

template<typename T>
class MatrixNonOwned {

public:

	MatrixNonOwned(SizeType nrow,
	               SizeType ncol,
	               typename VectorConstOrNot<T>::Type& data,
	               SizeType offset)
	    : nrow_(nrow), ncol_(ncol), data_(&data), offset_(offset)
	{
		check();
	}

	explicit MatrixNonOwned(const PsimagLite::Matrix<typename RemoveConst<T>::Type>& m)
	    : nrow_(m.n_row()), ncol_(m.n_col()), data_(&(m.data_)), offset_(0)
	{
		check();
	}

	explicit MatrixNonOwned(PsimagLite::Matrix<typename RemoveConst<T>::Type>& m)
	    : nrow_(m.n_row()), ncol_(m.n_col()), data_(&(m.data_)), offset_(0)
	{
		check();
	}

	const T& operator()(SizeType i, SizeType j) const
	{
		assert(i < nrow_);
		assert(j < ncol_);
		assert(offset_ + i + j*nrow_ < data_->size());
		return (*data_)[offset_ + i + j*nrow_];
	}

	T& operator()(SizeType i, SizeType j)
	{
		assert(i < nrow_);
		assert(j < ncol_);
		assert(offset_ + i + j*nrow_ < data_->size());
		return (*data_)[offset_ + i + j*nrow_];
	}

	const typename VectorConstOrNot<T>::Type& getVector() const
	{
		assert(data_);
		return *data_;
	}

	typename VectorConstOrNot<T>::Type& getVector()
	{
		assert(data_);
		return *data_;
	}

private:

	void check() const
	{
		assert(data_);
		assert(offset_ + nrow_ - 1 + (ncol_ - 1)*nrow_ < data_->size());
	}

	MatrixNonOwned(const MatrixNonOwned&);

	MatrixNonOwned& operator=(const MatrixNonOwned&);

	SizeType nrow_;
	SizeType ncol_;
	typename VectorConstOrNot<T>::Type* data_;
	SizeType offset_;
};
}

#endif
