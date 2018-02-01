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

		data_.resize(partitions.size() -1, partitions_.size() - 1);
		data_.setTo(0);

		MatrixBlockType value(size_, size_);
		value.setTo(0.0);

		VectorSizeType indexToPart(size_, 0);
		fillIndexToPart(indexToPart);
		PsimagLite::Matrix<SizeType> icount(n, n);
		icount.setTo(0);

		for (SizeType row = 0; row < size_; ++row) {
			SizeType kStart = sparse.getRowPtr(row);
			SizeType kEnd = sparse.getRowPtr(row + 1);
			SizeType ipatch = indexToPart[row];
			for (SizeType k = kStart; k < kEnd; ++k) {
				SizeType col = sparse.getCol(k);
				SizeType jpatch = indexToPart[col];
				icount(ipatch, jpatch)++;
			}
		}

		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			SizeType rows = partitions_[ipatch + 1] - partitions_[ipatch];
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				SizeType cols = partitions_[jpatch + 1] - partitions_[jpatch];
				if (icount(ipatch, jpatch) == 0) continue;
				data_(ipatch, jpatch) = new MatrixBlockType(rows, cols);
			}
		}

		for (SizeType row = 0; row < size_; ++row) {
			SizeType kStart = sparse.getRowPtr(row);
			SizeType kEnd = sparse.getRowPtr(row + 1);
			SizeType ipatch = indexToPart[row];
			for (SizeType k = kStart; k < kEnd; ++k) {
				SizeType col = sparse.getCol(k);
				SizeType jpatch = indexToPart[col];
				MatrixBlockType& m = *(data_(ipatch, jpatch));
				m(row - partitions_[ipatch], col - partitions_[jpatch]) = sparse.getValue(k);
			}
		}
	}

	~BlockOffDiagMatrix()
	{
		assert(partitions_.size() > 0);
		SizeType n = partitions_.size() - 1;
		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				MatrixBlockType* m = data_(ipatch, jpatch);
				delete m;
				data_(ipatch, jpatch) = 0;
			}
		}
	}

	void toSparse(SparseMatrixType& sparse) const
	{
		assert(partitions_.size() > 0);
		SizeType n = partitions_.size() - 1;
		VectorSizeType nonzeroInThisRow(size_, 0);
		SizeType count = 0;
		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				const MatrixBlockType* mptr = data_(ipatch, jpatch);
				if (mptr == 0) continue;
				const MatrixBlockType& m = *mptr;
				for (SizeType r = 0; r < m.rows(); ++r) {
					SizeType row = r + partitions_[ipatch];
					nonzeroInThisRow[row] += m.cols();
					count += m.cols()*m.rows();
				}
			}
		}

		sparse.clear();
		sparse.resize(size_, size_, count);

		count = 0;
		for (SizeType row = 0; row < size_; ++row) {
			sparse.setRow(row, count);
			count += nonzeroInThisRow[row];
		}

		sparse.setRow(size_, count);

		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				const MatrixBlockType* mptr = data_(ipatch, jpatch);
				if (mptr == 0) continue;
				const MatrixBlockType& m = *mptr;
				for (SizeType r = 0; r < m.rows(); ++r) {
					SizeType ip = sparse.getRowPtr(r + partitions_[ipatch]);
					for (SizeType c = 0; c < m.cols(); ++c) {
						sparse.setValues(ip, m(r, c));
						sparse.setCol(ip, c + partitions_[jpatch]);
						++ip;
 					}
				}
			}
		}

		sparse.checkValidity();
	}

	void transform(const BlockDiagonalMatrixType& f)
	{
		assert(partitions_.size() > 0);
		SizeType n = partitions_.size() - 1;
		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				MatrixBlockType* mptr = data_(ipatch, jpatch);
				if (mptr == 0) continue;
				MatrixBlockType& m = *mptr;
				const MatrixBlockType& mRight = f(jpatch);
				const MatrixBlockType& mLeft = f(ipatch);

				if (mLeft.rows() == 0 || mRight.rows() == 0) {
					m.clear();
					continue;
				}

				assert(m.cols() == mRight.rows());
				assert(m.rows() == mLeft.rows());

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
	typename PsimagLite::Matrix<MatrixBlockType*> data_;
};
}
#endif // BLOCKOFFDIAGMATRIX_H
