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
	    : rows_(0), cols_(0)
	{
		storeOffsets(offsetRows_, lrs.left());
		storeOffsets(offsetCols_, lrs.right());
		rows_ = lrs.left().size();
		cols_ = lrs.right().size();

		SizeType ns = lrs.left().size();
		PackIndicesType packSuper(ns);
		SizeType offset = src.offset(srcIndex);
		SizeType targetQn = src.qn(srcIndex);
		GenIjPatchType g(lrs, targetQn);
		const VectorSizeType& patchesLeft = g(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = g(GenIjPatchType::RIGHT);
		SizeType npatches = patchesLeft.size();
		assert(npatches == patchesRight.size());

		data_.resize(npatches, 0);

		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType partL = patchesLeft[ipatch];
			SizeType partR = patchesRight[ipatch];
			SizeType offsetL = lrs.left().partition(partL);
			SizeType rtotal = lrs.left().partition(partL + 1) - offsetL;
			SizeType offsetR = lrs.right().partition(partR);
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

	void toVectorWithOffsets(VectorWithOffsetType& dest) const
	{
		err("BlockDiagWf::toVectorWithOffsets() not implemented yet\n");
	}

	void transform(char charLeft,
	               char charRight,
	               const BlockDiagonalMatrixType& tLeft,
	               const BlockDiagonalMatrixType& tRight)
	{
		err("BlockDiagWf::transform() not implemented yet\n");
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
		err("BlockDiagWf::storeOffsets() not implemented yet\n");
	}

	VectorSizeType offsetRows_;
	VectorSizeType offsetCols_;
	SizeType rows_;
	SizeType cols_;
	typename PsimagLite::Vector<MatrixType*>::Type data_;
};
}
#endif // BLOCKDIAGMATRIXWF_H
