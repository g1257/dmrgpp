#ifndef BLOCKDIAGMATRIXWF_H
#define BLOCKDIAGMATRIXWF_H
#include "CrsMatrix.h"
#include "BlockDiagonalMatrix.h"
#include "LAPACK.h"
#include "PackIndices.h"

namespace Dmrg {

template<typename GenIjPatchType, typename VectorWithOffsetType>
class BlockDiagWf {

	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType, SizeType> PairType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;
	typedef typename GenIjPatchType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisType BasisType;
	typedef PsimagLite::PackIndices PackIndicesType;

public:

	BlockDiagWf(const VectorWithOffsetType& src,
	            SizeType srcIndex,
	            const LeftRightSuperType& lrs)
	    : rows_(lrs.left().size()),
	      cols_(lrs.right().size()),
	      offsetRows_(lrs.left().partition()),
	      offsetCols_(lrs.right().partition()),
	      genIjPatch_(lrs, src.qn(srcIndex))
	{
		storeOffsets(offsetRows_, lrs.left());
		storeOffsets(offsetCols_, lrs.right());

		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		SizeType offset = src.offset(srcIndex);
		const VectorSizeType& patchesLeft = genIjPatch_(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = genIjPatch_(GenIjPatchType::RIGHT);
		SizeType npatches = patchesLeft.size();
		assert(npatches == patchesRight.size());

		data_.resize(npatches, 0);

		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType partL = patchesLeft[ipatch];
			SizeType partR = patchesRight[ipatch];
			SizeType offsetL = offsetRows_[partL];
			SizeType rtotal = offsetRows_[partL + 1] - offsetL;
			SizeType offsetR = offsetCols_[partR];
			SizeType ctotal = offsetCols_[partR + 1] - offsetR;
			data_[ipatch] = new MatrixType(rtotal, ctotal);
			MatrixType& m = *(data_[ipatch]);
			for (SizeType r = 0; r < rtotal; ++r) {
				SizeType row = r + offsetL;
				for (SizeType c = 0; c < ctotal; ++c) {
					SizeType col = c + offsetR;
					SizeType ind = packSuper.pack(row,
					                              col,
					                              lrs.super().permutationInverse());
					assert(ind >= offset);
					m(r, c) = src.fastAccess(srcIndex, ind - offset);
				}
			}
		}
	}

	~BlockDiagWf()
	{
		SizeType npatches = data_.size();
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			delete data_[ipatch];
			data_[ipatch] = 0;
		}
	}

	void transform(char charLeft,
	               char charRight,
	               const BlockDiagonalMatrixType& tLeft,
	               const BlockDiagonalMatrixType& tRight)
	{
		const VectorSizeType& patchesLeft = genIjPatch_(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = genIjPatch_(GenIjPatchType::RIGHT);

		SizeType npatches = data_.size();
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			MatrixType* mptr = data_[ipatch];

			if (mptr == 0) continue;

			MatrixType& m = *mptr;
			const MatrixType& mRight = tRight(patchesRight[ipatch]);
			const MatrixType& mLeft = tLeft(patchesLeft[ipatch]);
			MatrixType tmp(m.rows(), (charRight == 'N') ? mRight.cols() : mRight.rows());

			if (mRight.cols() == 0 && mRight.rows() == 0) {
				m.clear();
				continue;
			}

			if (mLeft.cols() == 0 && mLeft.rows() == 0) {
				m.clear();
				continue;
			}

			// tmp = data_[ii] * op(mRight, charRight);

			// columns of the matrix op(A) ==  rows of the matrix op(B)
			assert(m.cols() == (charRight == 'N') ? mRight.rows() : mRight.cols());

			psimag::BLAS::GEMM('N',
			                   charRight,
			                   m.rows(),
			                   tmp.cols(),
			                   m.cols(),
			                   1.0,
			                   &(m(0,0)),
			                   m.rows(),
			                   &(mRight(0,0)),
			                   mRight.rows(),
			                   0.0,
			                   &(tmp(0,0)),
			                   tmp.rows());

			// data_[ii] = op(mLeft, charLeft) * tmp;

			// columns of the matrix op(A) ==  rows of the matrix op(B)
			assert((charLeft == 'N') ? mLeft.cols() : mLeft.rows() == tmp.rows());

			m.clear();
			m.resize((charLeft == 'N') ? mLeft.rows() : mLeft.cols(), tmp.cols());
			psimag::BLAS::GEMM(charLeft,
			                   'N',
			                   m.rows(),
			                   m.cols(),
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

		offsetRows_ = (charLeft == 'N')  ? tLeft.offsetsRows() : tLeft.offsetsCols();
		offsetCols_ = (charRight == 'N') ? tRight.offsetsCols() : tRight.offsetsRows();
		rows_ = (charLeft == 'N')  ? tLeft.rows() : tLeft.cols();
		cols_ = (charRight == 'N') ? tRight.cols() : tRight.rows();
	}

	void toVectorWithOffsets(VectorWithOffsetType& dest,
	                         SizeType destIndex,
	                         const LeftRightSuperType& lrs) const
	{
		const VectorSizeType& patchesLeft = genIjPatch_(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = genIjPatch_(GenIjPatchType::RIGHT);
		SizeType npatches = patchesLeft.size();

		assert(npatches == patchesRight.size());

		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		SizeType offset = dest.offset(destIndex);
		SizeType max = dest.offset(destIndex + 1);
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType partL = patchesLeft[ipatch];
			SizeType partR = patchesRight[ipatch];
			SizeType offsetL = offsetRows_[partL];
			SizeType offsetR = offsetCols_[partR];
			const MatrixType* mptr = data_[ipatch];

			if (mptr == 0) continue;

			const MatrixType& m = *mptr;
			for (SizeType r = 0; r < m.rows(); ++r) {
				SizeType row = r + offsetL;
				for (SizeType c = 0; c < m.cols(); ++c) {
					SizeType col = c + offsetR;
					SizeType ind = packSuper.pack(row,
					                              col,
					                              lrs.super().permutationInverse());
					if (ind < offset || ind >= max) continue;
					dest.fastAccess(destIndex, ind - offset) = m(r, c);
				}
			}
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

private:

	static void storeOffsets(VectorSizeType& offsets, const BasisType& basis)
	{
		SizeType n = basis.partition();
		if (n != offsets.size())
			err("BlockDiagWf::storeOffsets() offsets size\n");

		for (SizeType i = 0; i < n; ++i)
			offsets[i] = basis.partition(i);
	}

	SizeType rows_;
	SizeType cols_;
	VectorSizeType offsetRows_;
	VectorSizeType offsetCols_;
	GenIjPatchType genIjPatch_;
	typename PsimagLite::Vector<MatrixType*>::Type data_;
};
}
#endif // BLOCKDIAGMATRIXWF_H
