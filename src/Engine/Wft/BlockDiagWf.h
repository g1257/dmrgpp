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
	            SizeType srcIndex,
	            const LeftRightSuperType& lrs)
	    : lrs_(lrs),
	      rows_(lrs.left().size()),
	      cols_(lrs.right().size()),
	      offsetRows_(lrs.left().partition()),
	      offsetCols_(lrs.right().partition())
	{
		storeOffsets(offsetRows_, lrs.left());
		storeOffsets(offsetCols_, lrs.right());

		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		SizeType offset = src.offset(srcIndex);
		GenIjPatchType genIjPatch(lrs, src.qn(srcIndex));
		const VectorSizeType& patchesLeft = genIjPatch(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = genIjPatch(GenIjPatchType::RIGHT);
		SizeType npatches = patchesLeft.size();
		assert(npatches == patchesRight.size());

		data_.resize(npatches, 0);
		qns_.resize(npatches);
		patches_.resize(npatches);

		ComplexOrRealType sum = 0.0;
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType partL = patchesLeft[ipatch];
			SizeType partR = patchesRight[ipatch];
			SizeType offsetL = offsetRows_[partL];
			SizeType qnL = lrs.left().qn(offsetL);
			SizeType rtotal = offsetRows_[partL + 1] - offsetL;
			SizeType offsetR = offsetCols_[partR];
			SizeType qnR = lrs.right().qn(offsetR);
			qns_[ipatch] = PairType(qnL, qnR);
			patches_[ipatch] = PairType(partL, partR);
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
					sum += PsimagLite::conj(m(r, c))*m(r, c);
				}
			}
		}

		std::cout<<"sum -----------> "<<sum<<"\n";
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

		patchConvert(patchConvertLeft, false, tLeft);
		patchConvert(patchConvertRight, true, tRight);

		SizeType npatches = data_.size();
		ComplexOrRealType sum = 0.0;
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

			const MatrixType& mRight = tRight(partR);
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

			sum += normOfMatrix(m);

		}

		offsetRows_ = (charLeft == 'N')  ? tLeft.offsetsRows() : tLeft.offsetsCols();
		offsetCols_ = (charRight == 'N') ? tRight.offsetsCols() : tRight.offsetsRows();
		rows_ = (charLeft == 'N')  ? tLeft.rows() : tLeft.cols();
		cols_ = (charRight == 'N') ? tRight.cols() : tRight.rows();
		std::cout<<"sum transform "<<sum<<"\n";
	}

	void toVectorWithOffsets(VectorWithOffsetType& dest,
	                         SizeType destIndex,
	                         const LeftRightSuperType& lrs) const
	{
		VectorSizeType patchesLeft(qns_.size(), 0);
		VectorSizeType patchesRight(qns_.size(), 0);

		getPatchesLocations(patchesLeft, lrs.left(), true);
		getPatchesLocations(patchesRight, lrs.right(), false);

		SizeType npatches = patchesLeft.size();

		assert(npatches == patchesRight.size());

		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		//SizeType offset = lrs.super().partition(destIndex);
		//SizeType max = lrs.super().partition(destIndex + 1);
		ComplexOrRealType sum = 0.0;
		VectorType v(lrs.super().size(), 0.0);
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType partL = patchesLeft[ipatch];
			if (partL >= lrs.left().partition()) continue;
			SizeType partR = patchesRight[ipatch];
			if (partR >= lrs.right().partition()) continue;
			SizeType offsetL = lrs.left().partition(partL);
			SizeType offsetR = lrs.right().partition(partR);
			const MatrixType* mptr = data_[ipatch];

			if (mptr == 0) continue;

			const MatrixType& m = *mptr;

			for (SizeType r = 0; r < m.rows(); ++r) {
				SizeType row = r + offsetL;
				if (row >= lrs.left().size())
					break;
				for (SizeType c = 0; c < m.cols(); ++c) {
					SizeType col = c + offsetR;
					if (col >= lrs.right().size())
						break;
					SizeType ind = packSuper.pack(row,
					                              col,
					                              lrs.super().permutationInverse());
					v[ind] = m(r, c);
					sum += PsimagLite::conj(v[ind])*v[ind];
					//if (ind < offset || ind >= max) continue;
					//assert(ind >= offset && ind < max);
					//dest.fastAccess(destIndex, ind - offset) = m(r, c);
				}
			}
		}

		std::cout<<"sum = "<<sum<<"\n";
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

	void getPatchesLocations(VectorSizeType& patches,
	                         const BasisType& basis,
	                         bool isFirst) const
	{
		SizeType npatches = qns_.size();
		assert(patches.size() == npatches);
		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType qn = (isFirst) ? qns_[ipatch].first : qns_[ipatch].second;
			patches[ipatch] = findPatch(basis, qn);
		}
	}

	static SizeType findPatch(const BasisType& basis, SizeType q)
	{
		SizeType n = basis.partition();
		if (n == 0)
			err("BlockDiagWf::findPatch(): no partitions in basis\n");
		--n;
		for (SizeType i = 0; i < n; ++i) {
			SizeType offset = basis.partition(i);
			if (static_cast<SizeType>(basis.qn(offset)) == q) return i;
		}

		return n + 1;
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

	const LeftRightSuperType& lrs_;
	SizeType rows_;
	SizeType cols_;
	VectorSizeType offsetRows_;
	VectorSizeType offsetCols_;
	VectorPairType qns_;
	VectorPairType patches_;
	typename PsimagLite::Vector<MatrixType*>::Type data_;
};
}
#endif // BLOCKDIAGMATRIXWF_H
