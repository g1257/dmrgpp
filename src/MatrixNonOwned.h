#ifndef DMRG_MATRIX_NON_OWNED_H
#define DMRG_MATRIX_NON_OWNED_H
#include "Matrix.h"

namespace PsimagLite {

template<typename T>
class MatrixNonOwned {

public:

	MatrixNonOwned(SizeType nrow, SizeType ncol, T* data)
	    : nrow_(nrow), ncol_(ncol), data_(data)
	{}

	explicit MatrixNonOwned(const PsimagLite::Matrix<typename RemoveConst<T>::Type>& m)
	    : nrow_(m.n_row()), ncol_(m.n_col()), data_(&(m.data_[0]))
	{}

	explicit MatrixNonOwned(PsimagLite::Matrix<typename RemoveConst<T>::Type>& m)
	    : nrow_(m.n_row()), ncol_(m.n_col()), data_(&(m.data_[0]))
	{}

	const T& operator()(SizeType i, SizeType j) const
	{
//		assert(i + j*nrow_ < data_.size());
		return data_[i + j*nrow_];
	}

	T& operator()(SizeType i, SizeType j)
	{
//		assert(i + j*nrow_ < data_.size());
		return data_[i + j*nrow_];
	}

private:

	MatrixNonOwned(const MatrixNonOwned&);

	MatrixNonOwned& operator=(const MatrixNonOwned&);

	SizeType nrow_;
	SizeType ncol_;
	T* data_;
};
}

#endif
