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
	typedef PsimagLite::Vector<PairType>::Type VectorPairType;
	typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;
	typedef typename GenIjPatchType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisType BasisType;
	typedef PsimagLite::PackIndices PackIndicesType;

public:

	BlockDiagWf(const VectorWithOffsetType& src,
	            SizeType iSrc,
	            const LeftRightSuperType& lrs)
	    : lrs_(lrs),
	      rows_(lrs.left().size()),
	      cols_(lrs.right().size())
	{
		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		SizeType srcIndex = src.sector(iSrc);
		SizeType offset = src.offset(srcIndex);
		GenIjPatchType genIjPatch(lrs, src.qn(iSrc));
		const VectorSizeType& patchesLeft = genIjPatch(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = genIjPatch(GenIjPatchType::RIGHT);
		SizeType npatches = patchesLeft.size();
		assert(npatches == patchesRight.size());

		data_.resize(npatches, 0);
		patches_.resize(npatches);

		//ComplexOrRealType sum = 0.0;
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType partL = patchesLeft[ipatch];
			SizeType partR = patchesRight[ipatch];
			SizeType offsetL = lrs.left().partition(partL);
			SizeType rtotal = lrs.left().partition(partL + 1) - offsetL;
			SizeType offsetR = lrs.right().partition(partR);
			patches_[ipatch] = PairType(partL, partR);
			SizeType ctotal = lrs.right().partition(partR + 1) - offsetR;
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
					//sum += PsimagLite::conj(m(r, c))*m(r, c);
				}
			}
		}

		// std::cout<<"sum -----------> "<<sum<<"\n";
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
		VectorSizeType patchConvertLeft(tLeft.blocks(), 0);
		VectorSizeType patchConvertRight(tRight.blocks(), 0);

		patchConvert(patchConvertLeft, (charLeft == 'N'), tLeft);
		patchConvert(patchConvertRight, (charRight != 'N'), tRight);

		SizeType npatches = data_.size();
		//ComplexOrRealType sum = 0.0;
		SizeType rowsum = 0;
		SizeType colsum = 0;
		offsetCols_.clear();
		offsetRows_.clear();
		offsetCols_.resize(npatches, 0);
		offsetRows_.resize(npatches, 0);
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			MatrixType* mptr = data_[ipatch];

			if (mptr == 0) continue;

			MatrixType& m = *mptr;
			SizeType partL = patchConvertLeft[patches_[ipatch].first];
			SizeType partR = patchConvertRight[patches_[ipatch].second];

			if (partL >= tLeft.blocks()) {
				m.clear();
				continue;
			}


			if (partR >= tRight.blocks()) {
				m.clear();
				continue;
			}

			offsetRows_[ipatch] = (charLeft == 'N') ? tLeft.offsetsRows(partL) :
			                                          tLeft.offsetsCols(partL);
			offsetCols_[ipatch] = (charRight == 'N') ? tRight.offsetsCols(partR) :
			                                           tRight.offsetsRows(partR);

			const MatrixType& mRight = getRightMatrix(tRight(partR), charRight);
			const MatrixType& mLeft = tLeft(partL);
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
			if (charRight == 'N')
				assert(m.cols() == mRight.rows());
			else
				assert(m.cols() == mRight.cols());

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

			m.clear();

			// columns of the matrix op(A) ==  rows of the matrix op(B)
			if (charLeft == 'N') {
				m.resize(mLeft.rows(), tmp.cols());
				assert(mLeft.cols() == tmp.rows());
			} else {
				m.resize(mLeft.cols(), tmp.cols());
				assert(mLeft.rows() == tmp.rows());
			}

			rowsum += m.rows();
			colsum += m.cols();

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

			//sum += normOfMatrix(m);

		}

		rows_ = (charLeft == 'N')  ? tLeft.rows() : tLeft.cols();
		cols_ = (charRight == 'N') ? tRight.cols() : tRight.rows();
		//std::cout<<"sum transform "<<sum<<" rowsum="<<rowsum<<" colsum="<<colsum<<"\n";
	}

	void toVectorWithOffsets(VectorWithOffsetType& dest,
	                         SizeType iNew,
	                         const LeftRightSuperType& lrs,
	                         const VectorSizeType& nk,
	                         typename ProgramGlobals::DirectionEnum dir) const
	{
		SizeType destIndex = dest.sector(iNew);
		if (dir == ProgramGlobals::EXPAND_SYSTEM)
			return toVectorExpandSys(dest, destIndex, lrs, nk);

		toVectorExpandEnviron(dest, destIndex, lrs, nk);
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

	void toVectorExpandSys(VectorWithOffsetType& dest,
	                       SizeType destIndex,
	                       const LeftRightSuperType& lrs,
	                       const VectorSizeType& nk) const
	{
		assert(nk.size() > 0);
		SizeType hilbert = nk[0];

		SizeType npatches = data_.size();
		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		PackIndicesType packLeft(lrs.left().size()/hilbert);
		PackIndicesType packRight(hilbert);
		SizeType offset = lrs.super().partition(destIndex);
		//ComplexOrRealType sum = 0.0;
		//ComplexOrRealType sumBad = 0.0;
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			const MatrixType* mptr = data_[ipatch];

			if (mptr == 0) continue;

			const MatrixType& m = *mptr;
			SizeType offsetL = offsetRows_[ipatch];
			SizeType offsetR = offsetCols_[ipatch];

			for (SizeType r = 0; r < m.rows(); ++r) {
				SizeType row = r + offsetL;
				for (SizeType c = 0; c < m.cols(); ++c) {
					SizeType col = c + offsetR;
					SizeType k = 0;
					SizeType rind = 0;
					packRight.unpack(k, rind, lrs_.right().permutation(col));
					assert(k < hilbert);
					SizeType lind = packLeft.pack(row, k, lrs.left().permutationInverse());

					SizeType ind = packSuper.pack(lind,
					                              rind,
					                              lrs.super().permutationInverse());
					const ComplexOrRealType& value = m(r, c);
					//sum += PsimagLite::conj(value)*value;
					//if (ind < offset || ind >= lrs.super().partition(destIndex + 1))
					//	sumBad += PsimagLite::conj(value)*value;
					assert(ind >= offset && ind < lrs.super().partition(destIndex + 1));
					dest.fastAccess(destIndex, ind - offset) = value;
				}
			}
		}

		//std::cout<<"sum = "<<sum<<" sumBad= "<<sumBad<<"\n";
	}

	void toVectorExpandEnviron(VectorWithOffsetType& dest,
	                           SizeType destIndex,
	                           const LeftRightSuperType& lrs,
	                           const VectorSizeType& nk) const
	{
		assert(nk.size() > 0);
		SizeType hilbert = nk[0];

		SizeType npatches = data_.size();
		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		PackIndicesType packLeft(lrs_.left().permutationInverse().size()/hilbert);
		PackIndicesType packRight(hilbert);
		SizeType offset = lrs.super().partition(destIndex);
		//ComplexOrRealType sum = 0.0;
		//ComplexOrRealType sumBad = 0.0;
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			const MatrixType* mptr = data_[ipatch];

			if (mptr == 0) continue;

			const MatrixType& m = *mptr;
			SizeType offsetL = offsetRows_[ipatch];
			SizeType offsetR = offsetCols_[ipatch];

			for (SizeType r = 0; r < m.rows(); ++r) {
				SizeType row = r + offsetL;
				for (SizeType c = 0; c < m.cols(); ++c) {
					SizeType col = c + offsetR;
					SizeType k = 0;
					SizeType lind = 0;
					packLeft.unpack(lind, k, lrs_.left().permutation(row));
					assert(k < hilbert);
					SizeType rind = packRight.pack(k, col, lrs.right().permutationInverse());

					SizeType ind = packSuper.pack(lind,
					                              rind,
					                              lrs.super().permutationInverse());
					const ComplexOrRealType& value = m(r, c);
					//sum += PsimagLite::conj(value)*value;
					//if (ind < offset || ind >= lrs.super().partition(destIndex + 1))
					//	sumBad += PsimagLite::conj(value)*value;
					assert(ind >= offset && ind < lrs.super().partition(destIndex + 1));
					dest.fastAccess(destIndex, ind - offset) = value;
				}
			}
		}

		//std::cout<<"sum = "<<sum<<" sumBad= "<<sumBad<<"\n";
	}

	static void patchConvert(VectorSizeType& v,
	                         bool isNeeded,
	                         const BlockDiagonalMatrixType& b)
	{
		SizeType n = b.blocks();
		assert(v.size() == n);
		SizeType c = 0;
		for (SizeType i = 0; i < n; ++i) {
			const MatrixType& m = b(i);
			if (m.rows() == 0 && isNeeded) {
				if (m.cols() != 0)
					err("patchConvert error\n");
				continue;
			}

			v[c++] = i;
		}
	}

	static ComplexOrRealType normOfMatrix(const MatrixType& m)
	{
		ComplexOrRealType sum = 0.0;
		for (SizeType j = 0; j < m.cols(); ++j)
			for (SizeType i = 0; i < m.rows(); ++i)
				sum += PsimagLite::conj(m(i, j))*m(i, j);

		return sum;
	}

	const MatrixType& getRightMatrix(const MatrixType& m, char c)
	{
		if (c != 'N') return m;

		return getRightMatrixT(m);
	}

	const MatrixType& getRightMatrixT(const PsimagLite::Matrix<std::complex<double> >& m)
	{
		storage_.clear();
		SizeType rows = m.rows();
		SizeType cols = m.cols();
		storage_.resize(rows, cols);
		for (SizeType j = 0; j < cols; ++j)
			for (SizeType i = 0; i < rows; ++i)
				storage_(i, j) = PsimagLite::conj(m(i, j));

		return storage_;
	}

	const MatrixType& getRightMatrixT(const PsimagLite::Matrix<double>& m)
	{
		return m;
	}

	const LeftRightSuperType& lrs_;
	SizeType rows_;
	SizeType cols_;
	VectorSizeType offsetRows_;
	VectorSizeType offsetCols_;
	VectorPairType patches_;
	MatrixType storage_;
	typename PsimagLite::Vector<MatrixType*>::Type data_;
};
}
#endif // BLOCKDIAGMATRIXWF_H
