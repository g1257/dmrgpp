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
		if (sparse.rows() != sparse.cols())
			err("BlockOffDiagMatrix::ctor() expects square sparse matrix\n");

		if (partitions.size() == 0)
			err("BlockOffDiagMatrix::ctor() expects partitions.size() > 0\n");

		SizeType n = partitions.size() - 1;
		size_ = partitions_[n];

		MatrixBlockType value(size_, size_);
		value.setTo(0.0);

		VectorSizeType indexToPart(size_, 0);
		fillIndexToPart(indexToPart);
		VectorSizeType parts(n, 0);

		for (SizeType i = 0; i < n; ++i) {
			SizeType startRow = partitions_[i];
			SizeType rtotal = partitions_[i + 1] - startRow;
			SizeType count = 0;
			for (SizeType r = 0; r < rtotal; ++r) {
				SizeType kStart = sparse.getRowPtr(r + startRow);
				SizeType kEnd = sparse.getRowPtr(r + startRow + 1);
				for (SizeType k = kStart; k < kEnd; ++k) {
					SizeType col = sparse.getCol(k);
					SizeType j = indexToPart[col];
					parts[j] = 1;
					value(r, col) = sparse.getValue(k);
					++count;
				}
			}

			if (count == 0) continue;

			for (SizeType j = 0; j < n; ++j) {
				if (parts[j] == 0) continue;
				assert(j + 1 < partitions_.size());
				SizeType startCol = partitions_[j];
				SizeType ctotal = partitions_[j + 1] - startCol;
				MatrixBlockType* mptr = new MatrixBlockType(rtotal, ctotal);
				MatrixBlockType& m = *mptr;
				for (SizeType r = 0; r < rtotal; ++r)
					for (SizeType c = 0; c < ctotal; ++c)
						m(r, c) = value(r, c + startCol);

				data_.push_back(mptr);
				pairs_.push_back(PairType(i,j));
			}
		}
	}

	~BlockOffDiagMatrix()
	{
		SizeType n = data_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete data_[i];
			data_[i] = 0;
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
				assert(ii < data_.size());
				const MatrixBlockType& m = *(data_[ii]);
				assert(i < partitions_.size());
				SizeType start = partitions_[i];
				assert(row >= start);
				SizeType r = row - start;
				if (r >= m.cols()) continue;

				assert(j + 1 < partitions_.size());
				SizeType colStart = partitions_[j];
				SizeType ctotal = partitions_[j + 1] - colStart;

				if (ctotal > m.cols())
					ctotal = m.cols();

				for (SizeType c = 0; c < ctotal; ++c) {
					sparse.pushCol(c + colStart);
					sparse.pushValue(m(r, c));
					++count;
				}
			}
		}

		sparse.setRow(size_, count);
	}

	void transform(const BlockDiagonalMatrixType& f)
	{
		SizeType n = pairs_.size();
		assert(n == data_.size());
		for (SizeType ii = 0; ii < n; ++ii) {
			SizeType i = pairs_[ii].first;
			SizeType j = pairs_[ii].second;
			const MatrixBlockType& mRight = f(j);
			const MatrixBlockType& mLeft = f(i);
			MatrixBlockType& m = *(data_[ii]);
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

		partitions_ = f.offsetsCols();
		n = partitions_.size();
		assert(n > 0);
		--n;
		size_ = partitions_[n];
	}

	SizeType rows() const
	{
		return size_;
	}

	SizeType cols() const
	{
		return size_;
	}

private:

	void findIndicesForEachRow(VectorVectorSizeType& indices) const
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

	void fillIndexToPart(VectorSizeType& indexToPart) const
	{
		assert(size_ == indexToPart.size());
		SizeType n = partitions_.size();
		assert(n > 0);
		--n;
		for (SizeType i = 0; i < n; ++i) {
			SizeType start = partitions_[i];
			SizeType total = partitions_[i + 1] - start;
			for (SizeType r = 0; r < total; ++r)
				indexToPart[r + start] = i;
		}
	}

	VectorSizeType partitions_;
	SizeType size_;
	typename PsimagLite::Vector<MatrixBlockType*>::Type data_;
	typename PsimagLite::Vector<PairType>::Type pairs_;
};
}
#endif // BLOCKOFFDIAGMATRIX_H
