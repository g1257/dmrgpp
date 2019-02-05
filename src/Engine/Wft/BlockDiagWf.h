#ifndef BLOCKDIAGMATRIXWF_H
#define BLOCKDIAGMATRIXWF_H
#include "CrsMatrix.h"
#include "BlockDiagonalMatrix.h"
#include "LAPACK.h"
#include "PackIndices.h"

#include <iostream>
#include <iomanip>

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
	typedef typename PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;

	class ParallelBlockCtor {

	public:

		ParallelBlockCtor(const VectorSizeType& patcheLeft,
		                  const VectorSizeType& patchesRight,
		                  const LeftRightSuperType& lrs,
		                  const VectorWithOffsetType& src,
		                  SizeType iSrc,
		                  VectorPairType& patches,
		                  VectorMatrixType& data)
		    : patchesLeft_(patcheLeft),
		      patchesRight_(patchesRight),
		      lrs_(lrs),
		      packSuper_(lrs.left().size()),
		      src_(src),
		      srcIndex_(src.sector(iSrc)),
		      offset_(src.offset(srcIndex_)),
		      patches_(patches),
		      data_(data)
		{}

		SizeType tasks() const { return patchesLeft_.size(); }

		void doTask(SizeType ipatch, SizeType)
		{
			SizeType partL = patchesLeft_[ipatch];
			SizeType partR = patchesRight_[ipatch];
			SizeType offsetL = lrs_.left().partition(partL);
			SizeType rtotal = lrs_.left().partition(partL + 1) - offsetL;
			SizeType offsetR = lrs_.right().partition(partR);
			patches_[ipatch] = PairType(partL, partR);
			SizeType ctotal = lrs_.right().partition(partR + 1) - offsetR;
			data_[ipatch] = new MatrixType(rtotal, ctotal);
			MatrixType& m = *(data_[ipatch]);
			for (SizeType r = 0; r < rtotal; ++r) {
				SizeType row = r + offsetL;
				for (SizeType c = 0; c < ctotal; ++c) {
					SizeType col = c + offsetR;
					SizeType ind = packSuper_.pack(row,
					                               col,
					                               lrs_.super().permutationInverse());
					assert(ind >= offset_);
					m(r, c) = src_.fastAccess(srcIndex_, ind - offset_);
					//sum += PsimagLite::conj(m(r, c))*m(r, c);
				}
			}
		}

	private:

		const VectorSizeType& patchesLeft_;
		const VectorSizeType& patchesRight_;
		const LeftRightSuperType& lrs_;
		const PackIndicesType packSuper_;
		const VectorWithOffsetType& src_;
		SizeType srcIndex_;
		SizeType offset_;
		VectorPairType& patches_;
		VectorMatrixType& data_;
	};

	class ParallelBlockTransform {

	public:

		ParallelBlockTransform(const BlockDiagonalMatrixType& tLeft,
		                       const BlockDiagonalMatrixType& tRight,
		                       char charLeft,
		                       char charRight,
		                       SizeType threads,
		                       const VectorPairType& patches,
		                       VectorSizeType& offsetRows,
		                       VectorSizeType& offsetCols,
		                       VectorMatrixType& data)
		    : tLeft_(tLeft),
		      tRight_(tRight),
		      patchConvertLeft_(tLeft.blocks(), 0),
		      patchConvertRight_(tRight.blocks(), 0),
		      charLeft_(charLeft),
		      charRight_(charRight),
		      storage_(threads),
		      patches_(patches),
		      offsetRows_(offsetRows),
		      offsetCols_(offsetCols),
		      data_(data)
		{
			patchConvert(patchConvertLeft_, (charLeft == 'N'), tLeft);
			patchConvert(patchConvertRight_, (charRight != 'N'), tRight);
			SizeType npatches = data_.size();
			offsetCols_.clear();
			offsetRows_.clear();
			offsetCols_.resize(npatches, 0);
			offsetRows_.resize(npatches, 0);
		}

		SizeType tasks() const { return data_.size(); }

		void doTask(SizeType ipatch, SizeType threadNum)
		{
			const int idebug = 0;

			MatrixType* mptr = data_[ipatch];

			if (mptr == 0) return;

			MatrixType& m = *mptr;
			SizeType partL = patchConvertLeft_[patches_[ipatch].first];
			SizeType partR = patchConvertRight_[patches_[ipatch].second];

			if (partL >= tLeft_.blocks()) {
				m.clear();
				return;
			}


			if (partR >= tRight_.blocks()) {
				m.clear();
				return;
			}

			offsetRows_[ipatch] = (charLeft_ == 'N') ? tLeft_.offsetsRows(partL) :
			                                           tLeft_.offsetsCols(partL);
			offsetCols_[ipatch] = (charRight_ == 'N') ? tRight_.offsetsCols(partR) :
			                                            tRight_.offsetsRows(partR);

			const MatrixType& mRight = getRightMatrix(tRight_(partR), charRight_, threadNum);
			const MatrixType& mLeft = tLeft_(partL);

			if (mRight.cols() == 0 && mRight.rows() == 0) {
				m.clear();
				return;
			}

			if (mLeft.cols() == 0 && mLeft.rows() == 0) {
				m.clear();
				return;
			}


			//  ---------------------------------
			//  Here W_L = mLeft,  W_R = mRight
			//
			//  Need to evaluate
			//
			//  Ynew = opL(W_L) * Yold * opR( W_R )
			//
			//  but note that
			//  Ynew is over-written by Yold
			//
			//  where
			//
			//  opL(W_L) = W_L                  if charLeft_ == 'N'
			//  opL(W_L) = transpose(W_L)       if charLeft_ == 'T'
			//  opL(W_L) = conj(transpose(W_L)) if charLeft_ == 'C'
			//
			//  opL(W_R) = W_R                  if charRight_ == 'N'
			//  opL(W_R) = transpose(W_R)       if charRight_ == 'T'
			//  opL(W_R) = conj(transpose(W_R)) if charRight_ == 'C'
			//  ---------------------------------
			//
			const int nrow_W_L = mLeft.rows();
			const int ncol_W_L = mLeft.cols();
			const ComplexOrRealType *W_L = &(mLeft(0,0));
			const int ldW_L = nrow_W_L;

			const int nrow_W_R = mRight.rows();
			const int ncol_W_R = mRight.cols();
			const ComplexOrRealType *W_R = &(mRight(0,0));
			const int ldW_R = nrow_W_R;


			const bool has_work = (nrow_W_L >= 1) && (ncol_W_L >= 1) &&
			        (nrow_W_R >= 1) && (ncol_W_R >= 1);
			if (!has_work) {
				m.clear();
				return;
			};


			const int nrow_Yold = m.rows();
			const int ncol_Yold = m.cols();
			ComplexOrRealType *Yold = &(m(0,0));
			const int ldYold = nrow_Yold;

			const int nrow_Ynew = (charLeft_ == 'N') ? nrow_W_L : ncol_W_L;
			const int ncol_Ynew = (charRight_ == 'N') ? ncol_W_R : nrow_W_R;
			const int ldYnew = nrow_Ynew;

			//  ----------------------------------
			//  Note: Ynew is over-written by Yold
			//  so delay assign pointer value to Ynew
			//  ----------------------------------

			int nrow_Ytemp = 0;
			int ncol_Ytemp = 0;
			int ldYtemp = 0;
			// ---------------------------
			// Method 1:
			// (1) Ytemp = opL(W_L) * Yold
			// (2) Ynew = Ytemp * opR(W_R)
			//
			// Method 2:
			// (1) Ytemp = Yold * opR( W_R )
			// (2) Ynew = opL(W_L) * Ytemp
			// ---------------------------

			assert( (charLeft_ == 'N') || (charLeft_ == 'T') || (charLeft_ == 'C'));
			assert( (charRight_ == 'N') || (charRight_ == 'T') || (charRight_ == 'C'));

			nrow_Ytemp = (charLeft_ == 'N') ? nrow_W_L : ncol_W_L;
			ncol_Ytemp = ncol_Yold;
			const RealType flops_method_1 = 2.0 * nrow_W_L * ncol_W_L * ncol_Yold +
			        2.0 * nrow_Ytemp * ncol_Ytemp * ncol_Ynew;

			nrow_Ytemp = nrow_Yold;
			ncol_Ytemp = (charRight_ == 'N') ? ncol_W_R : nrow_W_R;
			const RealType flops_method_2 = 2.0 * nrow_Yold * ncol_Yold * ncol_Ytemp +
			        2.0 * nrow_W_L * ncol_W_L * ncol_Ytemp;

			const bool use_method_1 = (flops_method_1 <= flops_method_2);

			if (idebug >= 1) {
				const double speedup_ratio = (use_method_1) ?
				            flops_method_2 / flops_method_1 :
				            flops_method_1 / flops_method_2;

				std::cout << "BlockDiagWf.h:225: "
				          << " use_method_1=" << use_method_1
				          << std::scientific
				          << " flops_method_1=" << flops_method_1
				          << " flops_method_2=" << flops_method_2
				          << " speedup ratio=" << speedup_ratio
				          << std::defaultfloat
				          << "\n";
			};



			const ComplexOrRealType d_one = 1.0;
			const ComplexOrRealType d_zero = 0.0;

			if (use_method_1) {
				// ---------------------------
				// Method 1:
				// (1) Ytemp = opL(W_L) * Yold
				// (2) Ynew = Ytemp * opR(W_R)
				// ---------------------------
				nrow_Ytemp = (charLeft_ == 'N') ? nrow_W_L : ncol_W_L;
				ncol_Ytemp = ncol_Yold;
				MatrixType tmp(nrow_Ytemp,ncol_Ytemp);
				ComplexOrRealType *Ytemp = &(tmp(0,0));
				ldYtemp = nrow_Ytemp;

				// ---------------------------
				// (1) Ytemp = opL(W_L) * Yold
				// ---------------------------
				{
					const char transA = charLeft_;
					const char transB = 'N';
					const int mm = nrow_Ytemp;
					const int nn = ncol_Ytemp;
					const int kk = nrow_Yold;
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					psimag::BLAS::GEMM( transA, transB,
					                    mm, nn, kk,
					                    alpha,  W_L, ldW_L, Yold, ldYold,
					                    beta,   Ytemp, ldYtemp );
				}

				// ------------------------------
				// Note Ynew is over-written Yold
				// ------------------------------
				m.clear();
				m.resize( nrow_Ynew, ncol_Ynew);
				ComplexOrRealType  *Ynew = &(m(0,0));


				// ---------------------------
				// (2) Ynew = Ytemp * opR(W_R)
				// ---------------------------
				{
					const char transA = 'N';
					const char transB = charRight_;
					const int mm = nrow_Ynew;
					const int nn = ncol_Ynew;
					const int kk = ncol_Ytemp;
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					psimag::BLAS::GEMM( transA, transB,
					                    mm, nn, kk,
					                    alpha, Ytemp, ldYtemp, W_R, ldW_R,
					                    beta, Ynew, ldYnew );
				}
			}
			else {
				// ---------------------------
				// Method 2:
				// (1) Ytemp = Yold * opR( W_R )
				// (2) Ynew = opL(W_L) * Ytemp
				// ---------------------------

				nrow_Ytemp = nrow_Yold;
				ncol_Ytemp = (charRight_ == 'N') ? ncol_W_R : nrow_W_R;
				MatrixType tmp(nrow_Ytemp,ncol_Ytemp);
				ComplexOrRealType *Ytemp = &(tmp(0,0));
				ldYtemp = nrow_Ytemp;

				// ------------------------------
				// (1) Ytemp = Yold * opR( W_R )
				// ------------------------------
				{
					const char transA = 'N';
					const char transB = charRight_;
					const int mm = nrow_Ytemp;
					const int nn = ncol_Ytemp;
					const int kk = ncol_Yold;
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					psimag::BLAS::GEMM( transA, transB,
					                    mm, nn, kk,
					                    alpha, Yold, ldYold, W_R, ldW_R,
					                    beta,  Ytemp, ldYtemp );

				}

				// ------------------------------
				// Note Ynew over-written by Yold
				// ------------------------------
				m.clear();
				m.resize( nrow_Ynew, ncol_Ynew );
				ComplexOrRealType *Ynew = &(m(0,0));

				// ---------------------------
				// (2) Ynew = opL(W_L) * Ytemp
				// ---------------------------
				{
					const char transA = charLeft_;
					const char transB = 'N';
					const int mm = nrow_Ynew;
					const int nn = ncol_Ynew;
					const int kk = nrow_Ytemp;
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					psimag::BLAS::GEMM( transA, transB,
					                    mm, nn, kk,
					                    alpha, W_L, ldW_L, Ytemp, ldYtemp,
					                    beta, Ynew, ldYnew );
				}


			};

		}

	private:

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

		const MatrixType& getRightMatrix(const MatrixType& m, char c, SizeType threadNum)
		{
			if (c != 'N') return m;

			return getRightMatrixT(m, threadNum);
		}

		const MatrixType& getRightMatrixT(const PsimagLite::Matrix<std::complex<RealType> >& m,
		                                  SizeType threadNum)
		{
			storage_[threadNum].clear();
			SizeType rows = m.rows();
			SizeType cols = m.cols();
			storage_[threadNum].resize(rows, cols);
			for (SizeType j = 0; j < cols; ++j)
				for (SizeType i = 0; i < rows; ++i)
					storage_[threadNum](i, j) = PsimagLite::conj(m(i, j));

			return storage_[threadNum];
		}

		const MatrixType& getRightMatrixT(const PsimagLite::Matrix<RealType>& m, SizeType)
		{
			return m;
		}

		const BlockDiagonalMatrixType& tLeft_;
		const BlockDiagonalMatrixType& tRight_;
		VectorSizeType patchConvertLeft_;
		VectorSizeType patchConvertRight_;
		char charLeft_;
		char charRight_;
		typename PsimagLite::Vector<MatrixType>::Type storage_;
		const VectorPairType& patches_;
		VectorSizeType& offsetRows_;
		VectorSizeType& offsetCols_;
		VectorMatrixType& data_;
	};

public:

	BlockDiagWf(const VectorWithOffsetType& src,
	            SizeType iSrc,
	            const LeftRightSuperType& lrs)
	    : lrs_(lrs),
	      rows_(lrs.left().size()),
	      cols_(lrs.right().size())
	{
		GenIjPatchType genIjPatch(lrs, src.qn(iSrc));
		const VectorSizeType& patchesLeft = genIjPatch(GenIjPatchType::LEFT);
		const VectorSizeType& patchesRight = genIjPatch(GenIjPatchType::RIGHT);
		SizeType npatches = patchesLeft.size();
		assert(npatches == patchesRight.size());

		data_.resize(npatches, 0);
		patches_.resize(npatches);

		SizeType threads = std::min(npatches, PsimagLite::Concurrency::codeSectionParams.npthreads);
		typedef PsimagLite::Parallelizer<ParallelBlockCtor> ParallelizerType;
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType threadedCtor(codeSectionParams);

		ParallelBlockCtor helper(patchesLeft, patchesRight, lrs, src, iSrc, patches_, data_);

		threadedCtor.loopCreate(helper);
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
		SizeType npatches = data_.size();
		SizeType threads = std::min(npatches, PsimagLite::Concurrency::codeSectionParams.npthreads);
		typedef PsimagLite::Parallelizer<ParallelBlockTransform> ParallelizerType;
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType threadedTransform(codeSectionParams);

		ParallelBlockTransform helper(tLeft,
		                              tRight,
		                              charLeft,
		                              charRight,
		                              threads,
		                              patches_,
		                              offsetRows_,
		                              offsetCols_,
		                              data_);

		threadedTransform.loopCreate(helper);

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

	const LeftRightSuperType& lrs_;
	SizeType rows_;
	SizeType cols_;
	VectorSizeType offsetRows_;
	VectorSizeType offsetCols_;
	VectorPairType patches_;
	MatrixType storage_;
	VectorMatrixType data_;
};
}
#endif // BLOCKDIAGMATRIXWF_H
