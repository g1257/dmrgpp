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

	typedef MatrixBlockType value_type;

	template<typename SomeBasisType>
	BlockOffDiagMatrix(const SomeBasisType& rowBasis,
	                   const SomeBasisType& colBasis,
	                   const typename SomeBasisType::QnType& qtarget)
	{
		typedef typename SomeBasisType::VectorQnType::value_type QnType;
		SizeType rowPatches = rowBasis.partition();
		offsetRows_.resize(rowPatches);
		for (SizeType i = 0; i < rowPatches; ++i)
			offsetRows_[i] = rowBasis.partition(i);
		rows_ = offsetRows_[rowPatches - 1];

		SizeType colPatches = colBasis.partition();
		offsetCols_.resize(colPatches);
		for (SizeType i = 0; i < colPatches; ++i)
			offsetCols_[i] = colBasis.partition(i);
		cols_ = offsetCols_[colPatches - 1];
		data_.resize(rowPatches - 1, colPatches - 1, 0);

		for (SizeType i = 0; i < rowPatches - 1; ++i) {
			SizeType rows = offsetRows_[i + 1] - offsetRows_[i];
			QnType qrow = rowBasis.qnEx(i);
			for (SizeType j = 0; j < colPatches - 1; ++j) {
				QnType qcol = colBasis.qnEx(j);
				QnType q(qrow, qcol);
				if (q != qtarget) continue;
				SizeType cols = offsetCols_[j + 1] - offsetCols_[j];
				data_(i, j) = new MatrixBlockType(rows, cols);
			}
		}
	}

	BlockOffDiagMatrix(const SparseMatrixType& sparse,
	                   const VectorSizeType& partitions)
	    : offsetRows_(partitions), rows_(0), cols_(0)
	{
		if (sparse.rows() != sparse.cols())
			err("BlockOffDiagMatrix::ctor() expects square sparse matrix\n");

		if (partitions.size() == 0)
			err("BlockOffDiagMatrix::ctor() expects partitions.size() > 0\n");

		SizeType n = partitions.size() - 1;
		rows_ = cols_ = partitions[n];

		data_.resize(partitions.size() -1, partitions.size() - 1);
		data_.setTo(0);

		VectorSizeType indexToPart(rows_, 0);
		fillIndexToPart(indexToPart, partitions);
		PsimagLite::Matrix<SizeType> icount(n, n);
		icount.setTo(0);

		for (SizeType row = 0; row < rows_; ++row) {
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
			SizeType rows = partitions[ipatch + 1] - partitions[ipatch];
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				SizeType cols = partitions[jpatch + 1] - partitions[jpatch];
				if (icount(ipatch, jpatch) == 0) continue;
				data_(ipatch, jpatch) = new MatrixBlockType(rows, cols);
			}
		}

		for (SizeType row = 0; row < rows_; ++row) {
			SizeType kStart = sparse.getRowPtr(row);
			SizeType kEnd = sparse.getRowPtr(row + 1);
			SizeType ipatch = indexToPart[row];
			for (SizeType k = kStart; k < kEnd; ++k) {
				SizeType col = sparse.getCol(k);
				SizeType jpatch = indexToPart[col];
				MatrixBlockType& m = *(data_(ipatch, jpatch));
				m(row - partitions[ipatch], col - partitions[jpatch]) = sparse.getValue(k);
			}
		}
	}

	~BlockOffDiagMatrix()
	{
		assert(offsetRows_.size() > 0);
		SizeType nr = offsetRows_.size() - 1;
		SizeType nc = (offsetCols_.size() == 0) ? nr : offsetCols_.size();
		for (SizeType ipatch = 0; ipatch < nr; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < nc; ++jpatch) {
				MatrixBlockType* m = data_(ipatch, jpatch);
				delete m;
				data_(ipatch, jpatch) = 0;
			}
		}
	}

	ComplexOrRealType operator*(const BlockOffDiagMatrix& other) const
	{
		ComplexOrRealType sum = 0;
		SizeType nr = data_.rows();
		SizeType nc = data_.cols();
		if (nr != other.data_.rows() || nc != other.data_.cols())
			operatorFailed("*");
		for (SizeType ipatch = 0; ipatch < nr; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < nc; ++jpatch) {
				MatrixBlockType* m = data_(ipatch, jpatch);
				MatrixBlockType* o = other.data_(ipatch, jpatch);
				if (!m && !o) continue;
				const bool b1 = (m && !o);
				const bool b2 = (o && !m);
				if (b1 || b2)
					operatorFailed("*");
				sum += scalarProduct(*m, *o);
			}
		}

		return sum;
	}

	void randomize()
	{
		SizeType nr = data_.rows();
		SizeType nc = data_.cols();
		RealType sum = 0;
		for (SizeType ipatch = 0; ipatch < nr; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < nc; ++jpatch) {
				MatrixBlockType* m = data_(ipatch, jpatch);
				if (m == 0) continue;
				m->randomize();
				sum += PsimagLite::norm2(*m);
			}
		}

		assert(sum > 0);
		RealType factor = 1.0/sqrt(sum);
		for (SizeType ipatch = 0; ipatch < nr; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < nc; ++jpatch) {
				MatrixBlockType* m = data_(ipatch, jpatch);
				if (m == 0) continue;
				(*m) *= factor;
			}
		}
	}

	template<typename SomeBasisType>
	void fromMatrixColumn(const MatrixBlockType& src,
	                      SizeType col,
	                      const SomeBasisType& super,
	                      SizeType partitionIndex)
	{
		SizeType start = super.partition(partitionIndex);
		SizeType end = super.partition(partitionIndex + 1) - start;
		SizeType nr = data_.rows();
		SizeType nc = data_.cols();

		for (SizeType ipatch = 0; ipatch < nr; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < nc; ++jpatch) {
				MatrixBlockType* mptr = data_(ipatch, jpatch);
				if (mptr == 0) continue;
				MatrixBlockType& m = *mptr;
				SizeType rows = m.rows();
				SizeType cols = m.cols();
				for (SizeType i = 0; i < rows; ++i) {
					SizeType lindex = i + offsetRows_[ipatch];
					for (SizeType j = 0; j < cols; ++j) {
						SizeType rindex = j + offsetCols_[jpatch];
						SizeType sindex = super.permutationInverse(lindex + rindex*rows_);
						if (sindex < start || sindex >= end) continue;
						sindex -= start;
						m(i, j) = src(sindex, col);
					}
				}
			}
		}
	}

	void toSparse(SparseMatrixType& sparse) const
	{
		if (offsetCols_.size() != 0)
			err("BlockOffDiagMatrix::toSparse() only for square matrix\n");

		assert(offsetRows_.size() > 0);
		SizeType n = offsetRows_.size() - 1;
		VectorSizeType nonzeroInThisRow(rows_, 0);
		SizeType count = 0;
		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				const MatrixBlockType* mptr = data_(ipatch, jpatch);
				if (mptr == 0) continue;
				const MatrixBlockType& m = *mptr;
				count += m.cols()*m.rows();
				for (SizeType r = 0; r < m.rows(); ++r) {
					SizeType row = r + offsetRows_[ipatch];
					nonzeroInThisRow[row] += m.cols();
				}
			}
		}

		sparse.clear();
		sparse.resize(rows_, rows_, count);

		count = 0;
		for (SizeType row = 0; row < rows_; ++row) {
			sparse.setRow(row, count);
			count += nonzeroInThisRow[row];
		}

		sparse.setRow(rows_, count);

		for (SizeType ipatch = 0; ipatch < n; ++ipatch) {
			for (SizeType jpatch = 0; jpatch < n; ++jpatch) {
				const MatrixBlockType* mptr = data_(ipatch, jpatch);
				if (mptr == 0) continue;
				const MatrixBlockType& m = *mptr;
				for (SizeType r = 0; r < m.rows(); ++r) {
					SizeType ip = sparse.getRowPtr(r + offsetRows_[ipatch]);
					for (SizeType c = 0; c < m.cols(); ++c) {
						sparse.setValues(ip, m(r, c));
						sparse.setCol(ip, c + offsetRows_[jpatch]);
						++ip;
					}
				}
			}
		}

		sparse.checkValidity();
	}

	void transform(const BlockDiagonalMatrixType& f)
	{
		if (offsetCols_.size() != 0)
			err("BlockOffDiagMatrix::transform() only for square matrix\n");

		assert(offsetRows_.size() > 0);
		SizeType n = offsetRows_.size() - 1;
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

			offsetRows_ = f.offsetsCols();
			n = offsetRows_.size();
			assert(n > 0);
			--n;
			cols_ = rows_ = offsetRows_[n];
		}
	}

	SizeType rows() const
	{
		return rows_;
	}

	SizeType cols() const
	{
		return cols_;
	}

	const VectorSizeType& offsets(bool option) const
	{
		return (option) ? offsetRows_ : offsetCols_;
	}

	const RealType norm2() const
	{
		SizeType n = data_.rows();
		SizeType m = data_.cols();
		RealType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < m; ++j) {
				MatrixBlockType* mptr = data_(i, j);
				if (mptr == 0) continue;
				sum += PsimagLite::norm2(*mptr);
			}
		}

		return sum;
	}

	BlockOffDiagMatrix& operator*=(const ComplexOrRealType& value)
	{
		SizeType n = data_.rows();
		SizeType m = data_.cols();
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < m; ++j) {
				MatrixBlockType* mptr = data_(i, j);
				if (mptr == 0) continue;
				(*mptr) *= value;
			}
		}

		return *this;
	}

	BlockOffDiagMatrix& operator/=(const ComplexOrRealType& value)
	{
		return (*this) *= (1.0/value);
	}

	friend RealType norm(const BlockOffDiagMatrix& m)
	{
		return sqrt(m.norm2());
	}

private:

	static void fillIndexToPart(VectorSizeType& indexToPart,
	                            const VectorSizeType& partitions)
	{
		SizeType n = partitions.size();
		assert(n > 0);
		--n;
		assert(partitions[n] == indexToPart.size());
		for (SizeType i = 0; i < n; ++i) {
			SizeType start = partitions[i];
			SizeType total = partitions[i + 1] - start;
			for (SizeType r = 0; r < total; ++r)
				indexToPart[r + start] = i;
		}
	}

	void operatorFailed(PsimagLite::String what) const
	{
		err("BlockOffDiagMatrix::operator" + what + " failed\n");
	}

	BlockOffDiagMatrix(const BlockOffDiagMatrix&);

	BlockOffDiagMatrix& operator=(const BlockOffDiagMatrix&);

	VectorSizeType offsetRows_;
	VectorSizeType offsetCols_;
	SizeType rows_;
	SizeType cols_;
	typename PsimagLite::Matrix<MatrixBlockType*> data_;
};
}
#endif // BLOCKOFFDIAGMATRIX_H
