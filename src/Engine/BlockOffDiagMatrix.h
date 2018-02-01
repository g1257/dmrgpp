#ifndef BLOCKOFFDIAGMATRIX_H
#define BLOCKOFFDIAGMATRIX_H
#include "CrsMatrix.h"
#include "BlockDiagonalMatrix.h"
#include "LAPACK.h"

namespace Dmrg {

template<typename MatrixBlockType>
class BlockOffDiagMatrix {

	typedef typename MatrixBlockType::value_type ComplexOrRealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef BlockDiagonalMatrix<MatrixBlockType> BlockDiagonalMatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType, SizeType> PairType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

public:

	BlockOffDiagMatrix(const SparseMatrixType& sparse,
	                   const VectorSizeType& partitions)
	    : partitions_(partitions), size_(0)
	{
		if (partitions.size() == 0)
			err("BlockOffDiagMatrix has no partitions\n");

		SizeType n = partitions.size() - 1;
		size_ = partitions_[n];

		for (SizeType i = 0; i < n; ++i) {
			SizeType start = partitions_[i];
			SizeType end = partitions_[i + 1];
			MatrixBlockType* mptr = new MatrixBlockType(end - start, end - start);
			MatrixBlockType& m = *mptr;
			m.setTo(0.0);
			SizeType count = 0;
			for (SizeType row = start; row < end; ++row) {
				SizeType colStart = sparse.getRowPtr(row);
				SizeType colEnd = sparse.getRowPtr(row + 1);
				SizeType j = findPartition(colStart, partitions);
				// FIXME: Only one partition
				for (SizeType k = colStart; k < colEnd; ++k) {
					SizeType col = sparse.getCol(k);
					m(row - start, col - partitions[j]) = sparse.getValue(k);
					++count;
				}
			}

			if (count == 0) {
				delete mptr;
				mptr = 0;
				continue;
			}

			data_.push_back(mptr);
			pairs_.push_back(PairType(i,j));
		}
	}

	void toSparse(SparseMatrixType& sparse) const
	{
		SizeType count = 0;
		sparse.resize(size_, size_);
		VectorVectorSizeType indices(size_);
		findIndicesForEachRow(indices);
		for (SizeType row = 0; row < size_; ++row) {
			sparse.setRow(row, count);
			const VectorSizeType& indicesForThisRow = indices[row];
			SizeType n = indicesForThisRow.size();
			for (SizeType iii = 0; iii < n; ++iii) {
				SizeType ii = indicesForThisRow[iii];
				assert(ii < pairs_.size());
				SizeType i = pairs_[ii].first;
				SizeType j = pairs_[ii].second;
				assert(i < partitions_.size());
				SizeType start = partitions_[i];
				assert(j + 1 < partitions_.size());
				SizeType colStart = partitions_[j];
				SizeType colEnd = partitions_[j + 1];
				for (SizeType col = colStart; col < colEnd; ++col) {
					sparse.pushCol(col);
					assert(row >= start);
					assert(col >= colStart);
					sparse.pushValue(data_(row - start, col - colStart));
					++count;
				}
			}
		}

		sparse.setRow(size_, count);
	}

	void transform(BlockDiagonalMatrixType& f)
	{
		SizeType n = pairs_.size();
		assert(n == data_.size());
		for (SizeType ii = 0; ii < n; ++ii) {
			SizeType i = pairs_[ii].first;
			SizeType j = pairs_[ii].second;
			const MatrixBlockType& mRight = f(j);
			const MatrixBlockType& mLeft = f(i);
			MatrixBlockType& m = data_[ii];
			MatrixBlockType tmp(m.rows(), mRight.cols());
			// tmp = data_[ii] * mRight;
			psimag::BLAS::GEMM('N',
			                   'N',
			                   m.rows(),
			                   mRight.cols(),
			                   m.cols(),
			                   1.0,
			                   &(m(0,0)),
			                   m.rows(),
			                   &(mRight(0,0)),
			                   mRight.rows(),
			                   0.0,
			                   &(tmp(0,0)),
			                   tmp.rows());
			// data_[ii] = transposeConjugate(mLeft) * tmp;
			m.clear();
			m.resize(mLeft.cols(), mRight.cols());
			psimag::BLAS::GEMM('C',
			                   'N',
			                   mLeft.cols(),
			                   tmp.cols(),
			                   tmp.rows(),
			                   1.0,
			                   &(mLeft(0,0)),
			                   mLeft.rows(),
			                   &(tmp(0,0)),
			                   tmp.rows(),
			                   0.0,
			                   &(m(0,0)),
			                   m.rows());
		}
	}

	SizeType rows() const
	{
		return size_;
	}

	SizeType cols() const
	{
		return size_;
	}

	SizeType blocks() const
	{
		return data_.size();
	}

private:

	SizeType findIndicesForEachRow(VectorVectorSizeType& indices) const
	{
		SizeType n = pairs_.size();
		for (SizeType ii = 0; ii < n; ++ii) {
			SizeType i = pairs_[ii].first;
			assert(i + 1 < partitions_.size());
			SizeType start = partitions_[i];
			SizeType end = partitions_[i + 1];
			assert(start <= end);
			assert(end <= indices.size());
			for (SizeType row = start; row < end; ++row)
				indices[row].push_back(ii);
		}
	}

	SizeType size_;
	VectorSizeType partitions_;
	typename PsimagLite::Vector<MatrixBlockType*>::Type data_;
	typename PsimagLite::Vector<PairType>::Type pairs_;
};
}
#endif // BLOCKOFFDIAGMATRIX_H
