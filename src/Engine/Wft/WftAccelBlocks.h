#ifndef WFTACCELBLOCKS_H
#define WFTACCELBLOCKS_H
#include "Matrix.h"
#include "BLAS.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename WaveFunctionTransfBaseType>
class WftAccelBlocks {

	typedef typename WaveFunctionTransfBaseType::DmrgWaveStructType DmrgWaveStructType;
	typedef typename WaveFunctionTransfBaseType::WftOptionsType WftOptionsType;
	typedef typename WaveFunctionTransfBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename WaveFunctionTransfBaseType::VectorSizeType VectorSizeType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename WaveFunctionTransfBaseType::PackIndicesType PackIndicesType;

	class ParallelWftInBlocks {

	public:

		ParallelWftInBlocks(VectorMatrixType& result,
		                    const VectorMatrixType& psi,
		                    const MatrixType& ws,
		                    const MatrixType& we,
		                    SizeType volumeOfNk,
		                    SizeType sysOrEnv,
		                    SizeType threads)
		    : result_(result),
		      psi_(psi),
		      ws_(ws),
		      we_(we),
		      volumeOfNk_(volumeOfNk),
		      sysOrEnv_(sysOrEnv),
		      storage_(threads)
		{}

		SizeType tasks() const { return volumeOfNk_; }

		void doTask(SizeType kp, SizeType threadNum)
		{
			if (sysOrEnv_ == ProgramGlobals::SYSTEM)
				return doTaskSystem(kp);

			doTaskEnviron(kp, threadNum);
		}

	private:

		void doTaskEnviron(SizeType kp, SizeType threadNum)
		{
			SizeType ipsize = ws_.rows();
			SizeType i2psize = ws_.cols();
			SizeType jp2size = we_.rows();
			SizeType jpsize = we_.cols();
			MatrixType tmp(i2psize, jpsize);

			result_[kp].resize(ipsize, jpsize);
			result_[kp].setTo(0.0);
			tmp.setTo(0.0);

			const MatrixType& weModif = getWeModif(we_, threadNum);

			psimag::BLAS::GEMM('N',
			                   'N',
			                   i2psize,
			                   jpsize,
			                   jp2size,
			                   1.0,
			                   &((psi_[kp])(0,0)),
			                   i2psize,
			                   &(weModif(0,0)),
			                   jp2size,
			                   0.0,
			                   &(tmp(0,0)),
			                   i2psize);

			psimag::BLAS::GEMM('N',
			                   'N',
			                   ipsize,
			                   jpsize,
			                   i2psize,
			                   1.0,
			                   &(ws_(0,0)),
			                   ipsize,
			                   &(tmp(0,0)),
			                   i2psize,
			                   0.0,
			                   &((result_[kp])(0,0)),
			                   ipsize);
		}

		void doTaskSystem(SizeType kp)
		{
			SizeType ipSize = ws_.rows();
			SizeType isSize = ws_.cols();
			SizeType jenSize = we_.rows();
			SizeType jprSize = we_.cols();
			MatrixType tmp(ipSize, jenSize);

			result_[kp].resize(isSize, jenSize);
			result_[kp].setTo(0.0);
			tmp.setTo(0.0);

			psimag::BLAS::GEMM('N',
			                   'T',
			                   ipSize,
			                   jenSize,
			                   jprSize,
			                   1.0,
			                   &((psi_[kp])(0,0)),
			                   ipSize,
			                   &(we_(0,0)),
			                   jenSize,
			                   0.0,
			                   &(tmp(0,0)),
			                   ipSize);

			psimag::BLAS::GEMM('C',
			                   'N',
			                   isSize,
			                   jenSize,
			                   ipSize,
			                   1.0,
			                   &(ws_(0,0)),
			                   ipSize,
			                   &(tmp(0,0)),
			                   ipSize,
			                   0.0,
			                   &((result_[kp])(0,0)),
			                   isSize);
		}

		const MatrixType& getWeModif(const PsimagLite::Matrix<RealType>& m, SizeType)
		{
			return m;
		}

		const MatrixType& getWeModif(const PsimagLite::Matrix<std::complex<RealType> >& m,
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

		VectorMatrixType& result_;
		const VectorMatrixType& psi_;
		const MatrixType& ws_;
		const MatrixType& we_;
		SizeType volumeOfNk_;
		SizeType sysOrEnv_;
		VectorMatrixType storage_;
	};

public:

	WftAccelBlocks(const DmrgWaveStructType& dmrgWaveStruct,
	               const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct), wftOptions_(wftOptions)
	{}

	void environFromInfinite(VectorWithOffsetType& psiDest,
	                         SizeType i0,
	                         const VectorWithOffsetType& psiSrc,
	                         SizeType i0src,
	                         const LeftRightSuperType& lrs,
	                         const VectorSizeType& nk) const
	{
		if (lrs.left().block().size() < 2)
			err("Bounce!?\n");

		SizeType volumeOfNk = ProgramGlobals::volumeOf(nk);
		MatrixType ws;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SYSTEM).toDense(ws);

		MatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::ENVIRON).toDense(we);

		SizeType i2psize = ws.cols();
		SizeType jp2size = we.rows();

		VectorMatrixType psi(volumeOfNk);
		for (SizeType kp = 0; kp < volumeOfNk; ++kp) {
			psi[kp].resize(i2psize, jp2size);
			psi[kp].setTo(0.0);
		}

		environPreparePsi(psi, psiSrc, i0src, volumeOfNk);

		VectorMatrixType result(volumeOfNk);

		SizeType threads = std::min(volumeOfNk, PsimagLite::Concurrency::codeSectionParams.npthreads);
		typedef PsimagLite::Parallelizer<ParallelWftInBlocks> ParallelizerType;
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType threadedWft(codeSectionParams);

		ParallelWftInBlocks helperWft(result,
		                              psi,
		                              ws,
		                              we,
		                              volumeOfNk,
		                              ProgramGlobals::ENVIRON,
		                              threads);

		threadedWft.loopCreate(helperWft);

		environCopyOut(psiDest, i0, result, lrs, volumeOfNk);
	}

	void systemFromInfinite(VectorWithOffsetType& psiDest,
	                        SizeType i0,
	                        const VectorWithOffsetType& psiSrc,
	                        SizeType i0src,
	                        const LeftRightSuperType& lrs,
	                        const VectorSizeType& nk) const
	{
		if (lrs.right().block().size() < 2)
			err("Bounce!?\n");

		SizeType volumeOfNk = ProgramGlobals::volumeOf(nk);
		MatrixType ws;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SYSTEM).toDense(ws);

		MatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::ENVIRON).toDense(we);

		SizeType ipSize = ws.rows();
		SizeType jprSize = we.cols();

		VectorMatrixType psi(volumeOfNk);
		for (SizeType kp = 0; kp < volumeOfNk; ++kp) {
			psi[kp].resize(ipSize, jprSize);
			psi[kp].setTo(0.0);
		}

		systemPreparePsi(psi, psiSrc, i0src, volumeOfNk);

		VectorMatrixType result(volumeOfNk);

		SizeType threads = std::min(volumeOfNk, PsimagLite::Concurrency::codeSectionParams.npthreads);
		typedef PsimagLite::Parallelizer<ParallelWftInBlocks> ParallelizerType;
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType threadedWft(codeSectionParams);

		ParallelWftInBlocks helperWft(result,
		                              psi,
		                              ws,
		                              we,
		                              volumeOfNk,
		                              ProgramGlobals::SYSTEM,
		                              threads);

		threadedWft.loopCreate(helperWft);

		systemCopyOut(psiDest, i0, result, lrs, volumeOfNk);
	}

private:

	void environPreparePsi(VectorMatrixType& psi,
	                       const VectorWithOffsetType& psiSrc,
	                       SizeType i0src,
	                       SizeType volumeOfNk) const
	{
		SizeType total = psiSrc.effectiveSize(i0src);
		SizeType offset = psiSrc.offset(i0src);
		PackIndicesType packSuper(dmrgWaveStruct_.lrs().left().size());
		PackIndicesType packLeft(dmrgWaveStruct_.lrs().left().size()/volumeOfNk);

		for (SizeType x = 0; x < total; ++x) {
			SizeType alpha = 0;
			SizeType jp2 = 0;
			packSuper.unpack(alpha, jp2, dmrgWaveStruct_.lrs().super().permutation(x + offset));
			SizeType ip2 = 0;
			SizeType kp = 0;
			packLeft.unpack(ip2, kp, dmrgWaveStruct_.lrs().left().permutation(alpha));
			psi[kp](ip2, jp2) += psiSrc.fastAccess(i0src, x);
		}
	}

	void environCopyOut(VectorWithOffsetType& psiDest,
	                    SizeType i0,
	                    const VectorMatrixType& result,
	                    const LeftRightSuperType& lrs,
	                    SizeType volumeOfNk) const
	{
		SizeType nip = lrs.super().permutationInverse().size()/
		        lrs.right().permutationInverse().size();
		PackIndicesType pack1(nip);
		PackIndicesType pack2(volumeOfNk);
		SizeType total = psiDest.effectiveSize(i0);
		SizeType start = psiDest.offset(i0);

		for (SizeType x = 0; x < total; ++x) {
			SizeType ip = 0;
			SizeType beta = 0;
			pack1.unpack(ip, beta, lrs.super().permutation(x+start));
			SizeType kp = 0;
			SizeType jp = 0;
			pack2.unpack(kp, jp, lrs.right().permutation(beta));
			psiDest.fastAccess(i0, x) += result[kp](ip, jp);
		}
	}

	void systemPreparePsi(VectorMatrixType& psi,
	                      const VectorWithOffsetType& psiSrc,
	                      SizeType i0src,
	                      SizeType volumeOfNk) const
	{
		SizeType total = psiSrc.effectiveSize(i0src);
		SizeType offset = psiSrc.offset(i0src);
		PackIndicesType packSuper(dmrgWaveStruct_.lrs().left().size());
		PackIndicesType packRight(volumeOfNk);

		for (SizeType y = 0; y < total; ++y) {
			SizeType ip = 0;
			SizeType jp = 0;
			packSuper.unpack(ip, jp, dmrgWaveStruct_.lrs().super().permutation(y + offset));
			SizeType jpl = 0;
			SizeType jpr = 0;
			packRight.unpack(jpl, jpr, dmrgWaveStruct_.lrs().right().permutation(jp));
			psi[jpl](ip, jpr) = psiSrc.fastAccess(i0src, y);
		}
	}

	void systemCopyOut(VectorWithOffsetType& psiDest,
	                   SizeType i0,
	                   const VectorMatrixType& result,
	                   const LeftRightSuperType& lrs,
	                   SizeType volumeOfNk) const
	{
		SizeType nip = lrs.left().permutationInverse().size()/volumeOfNk;
		SizeType nalpha = lrs.left().permutationInverse().size();
		PackIndicesType pack1(nalpha);
		PackIndicesType pack2(nip);
		SizeType total = psiDest.effectiveSize(i0);
		SizeType start = psiDest.offset(i0);

		for (SizeType x = 0; x < total; ++x) {
			SizeType isn = 0;
			SizeType jen = 0;
			pack1.unpack(isn, jen, lrs.super().permutation(x+start));
			SizeType is = 0;
			SizeType jpl = 0;
			pack2.unpack(is, jpl, lrs.left().permutation(isn));
			psiDest.fastAccess(i0, x) += result[jpl](is, jen);
		}
	}

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
};
}
#endif // WFTACCELBLOCKS_H
