#ifndef WFTACCELBLOCKS_H
#define WFTACCELBLOCKS_H
#include "Matrix.h"
#include "BLAS.h"
#include "ProgramGlobals.h"

#include <iostream>
#include <iomanip>

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
		                    const ProgramGlobals::SysOrEnvEnum sysOrEnv,
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
			if (sysOrEnv_ == ProgramGlobals::SysOrEnvEnum::SYSTEM)
				return doTaskSystem(kp);

			doTaskEnviron(kp, threadNum);
		}

	private:

		void doTaskEnviron(SizeType kp, SizeType threadNum)
		{
#if (0)
			SizeType ipsize = ws_.rows();
			SizeType i2psize = ws_.cols();
			SizeType jp2size = we_.rows();
			SizeType jpsize = we_.cols();
#endif

			const int nrow_W_S = ws_.rows();
			const int ncol_W_S = ws_.cols();
			const int nrow_W_E = we_.rows();
			const int ncol_W_E = we_.cols();
			//
			// --------------------------------
			// Compute Ynew = W_S * Yold * W_E
			// --------------------------------
			const int nrow_Yold = ncol_W_S;
			const int ncol_Yold = nrow_W_E;

			const int nrow_Ynew = nrow_W_S;
			const int ncol_Ynew = ncol_W_E;

			result_[kp].resize(nrow_Ynew, ncol_Ynew);
			result_[kp].setTo(0.0);


			const MatrixType& weModif = getWeModif(we_, threadNum);

#if 0
			const int idebug = 0;
#endif

			int nrow_Ytemp = 0;
			int ncol_Ytemp = 0;


			// ---------------------------
			// Compute   Ynew = W_S * Yold * W_E
			// as
			// Method 1
			// (1) Ytemp = W_S * Yold
			// (2) Ynew = Ytemp * W_E
			// need   2 * nrow_W_S * ncol_W_S * ncol_Yold   flops +
			//        2 * nrow_Ytemp * ncol_Ytemp * ncol_W_E flops
			// or
			// Method 2:
			// (1) Ytemp = Yold * W_E
			// (2) Ynew = W_S * Ytemp
			// need 2 * nrow_Yold * ncol_Yold * ncol_W_E flops +
			//      2 * nrow_W_S * ncol_W_S * ncol_Ytemp flops
			// ---------------------------


			const ComplexOrRealType *Yold = &((psi_[kp])(0,0));
			const int ldYold = nrow_Yold;

			const ComplexOrRealType *W_E = &(weModif(0,0));
			const int ldW_E = nrow_W_E;

			const ComplexOrRealType *W_S = &(ws_(0,0));
			const int ldW_S = nrow_W_S;

			ComplexOrRealType *Ynew = &((result_[kp])(0,0));
			const int ldYnew = nrow_Ynew;

			// ----------------------
			// Method 1
			// (1) Ytemp = W_S * Yold
			// (2) Ynew = Ytemp * W_E
			// ----------------------
			nrow_Ytemp = nrow_W_S;
			ncol_Ytemp = ncol_Yold;
			const RealType flops_method_1 = 2.0 * nrow_W_S   * ncol_W_S   * ncol_Yold +
			        2.0 * nrow_Ytemp * ncol_Ytemp * ncol_W_E;

			// ------------------------
			// Method 2:
			// (1) Ytemp = Yold * W_E
			// (2) Ynew = W_S * Ytemp
			// ------------------------
			nrow_Ytemp = nrow_Yold;
			ncol_Ytemp = ncol_W_E;
			const RealType flops_method_2 =  2.0 * nrow_Yold * ncol_Yold * ncol_W_E +
			        2.0 * nrow_W_S  * ncol_W_S  * ncol_Ytemp ;


			const bool use_method_1 = (flops_method_1 <= flops_method_2);
#if 0
			if (idebug >= 1) {
				std::cout << "WftAccelBlocks.h:146: "
				          << " use_method_1=" << use_method_1
				          << std::scientific
				          << " flops_method_1=" << flops_method_1
				          << " flops_method_2=" << flops_method_2
				          << std::defaultfloat
				          << "\n";

			}
#endif

			const ComplexOrRealType d_one = 1.0;
			const ComplexOrRealType d_zero = 0.0;

			// ---------------------------
			// Compute   Ynew = W_S * (Yold * W_E)
			// ---------------------------
			//
			// Note Y = kron(A,B) * X = B * X * transpose(A)
			//
			// Ynew = kron( transpose(W_E), W_S) * Yold
			// ------------------------------------------

			// -----------------------------
			// Method 1: Ytemp = W_S * Yold
			// Method 2: Ytemp = Yold * W_E
			// -----------------------------
			nrow_Ytemp = (use_method_1) ? nrow_W_S : nrow_Yold;
			ncol_Ytemp = (use_method_1) ? ncol_Yold : ncol_W_E;

			MatrixType tmp(nrow_Ytemp, ncol_Ytemp);
			tmp.setTo(0.0);
			ComplexOrRealType *Ytemp = &(tmp(0,0));
			const int ldYtemp = nrow_Ytemp;

			if (use_method_1) {
				// ------------------
				// Method 1:
				// Ytemp = W_S * Yold
				// Ynew = Ytemp * W_E
				// ------------------

				// ------------------
				// Ytemp = W_S * Yold
				// ------------------
				{
					const int mm = nrow_Ytemp;
					const int nn = ncol_Ytemp;
					const int kk = nrow_Yold;

					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					assert( ncol_W_S == nrow_Yold );
					assert( nrow_Ytemp == nrow_W_S );
					assert( ncol_Ytemp == ncol_Yold );

					psimag::BLAS::GEMM( 'N', 'N',
					                    mm, nn, kk,
					                    alpha, W_S, ldW_S, Yold, ldYold,
					                    beta,  Ytemp, ldYtemp );
				}

				// -------------------
				// Ynew = Ytemp * W_E
				// -------------------
				{
					const int mm = nrow_Ynew;
					const int nn = ncol_Ynew;
					const int kk = ncol_Ytemp;

					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					assert( nrow_Ynew == nrow_Ytemp );
					assert( ncol_Ynew == ncol_W_E );
					assert( ncol_Ytemp == nrow_W_E );

					psimag::BLAS::GEMM( 'N', 'N',
					                    mm, nn, kk,
					                    alpha, Ytemp, ldYtemp, W_E, ldW_E,
					                    beta,  Ynew, ldYnew );


				}


			}
			else {
				// --------------------
				// Method 2:
				// Ytemp = Yold * W_E
				// Ynew =  W_S * Ytemp
				// --------------------


				//  ------------------
				//  Ytemp = Yold * W_E
				//  ------------------
				{
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					const int mm = nrow_Ytemp;
					const int nn = ncol_Ytemp;
					const int kk = ncol_Yold;

					assert( nrow_Ytemp == nrow_Yold );
					assert( ncol_Ytemp == ncol_W_E );
					assert( ncol_Yold == nrow_W_E );

					psimag::BLAS::GEMM('N', 'N',
					                   mm, nn, kk,
					                   alpha, Yold, ldYold, W_E, ldW_E,
					                   beta, Ytemp, ldYtemp );
				}

				// ------------------
				// Ynew = W_S * Ytemp
				// ------------------
				{
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					const int mm = nrow_Ynew;
					const int nn = ncol_Ynew;
					const int kk = ncol_W_S;

					assert( nrow_Ynew == nrow_W_S);
					assert( ncol_Ynew == ncol_Ytemp);
					assert( ncol_W_S == nrow_Ytemp );

					psimag::BLAS::GEMM('N', 'N',
					                   mm, nn, kk,
					                   alpha, W_S, ldW_S, Ytemp, ldYtemp,
					                   beta, Ynew, ldYnew );

				}
			};
		}

		void doTaskSystem(SizeType kp)
		{
#if (0)
			SizeType ipSize = ws_.rows();
			SizeType isSize = ws_.cols();
			SizeType jenSize = we_.rows();
			SizeType jprSize = we_.cols();
#endif

#if 0
			const int idebug = 0;
#endif

			const int nrow_W_S = ws_.rows();
			const int ncol_W_S = ws_.cols();
			const int nrow_W_E = we_.rows();
			const int ncol_W_E = we_.cols();



			// -------------------------------------------------
			// Note Y = kron(A,B) * X  =   B * X * transpose(A)
			//
			// Ynew = kron( W_E, conj(transpose(W_S)) * Yold
			// -------------------------------------------------
			const int nrow_Yold = nrow_W_S;
			const int ncol_Yold = ncol_W_E;
			const int nrow_Ynew = ncol_W_S;
			const int ncol_Ynew = nrow_W_E;

			int nrow_Ytemp = 0;
			int ncol_Ytemp = 0;


			result_[kp].resize(nrow_Ynew, ncol_Ynew);
			result_[kp].setTo(0.0);

			const ComplexOrRealType *Yold = &((psi_[kp])(0,0));
			const int ldYold = nrow_Yold;

			const ComplexOrRealType *W_E = &(we_(0,0));
			const int ldW_E = nrow_W_E;


			ComplexOrRealType  *Ynew = &((result_[kp])(0,0));
			const int ldYnew = nrow_Ynew;

			const ComplexOrRealType *W_S = &(ws_(0,0));
			const int ldW_S = nrow_W_S;


			const ComplexOrRealType d_one = 1.0;
			const ComplexOrRealType d_zero = 0.0;

			// ---------------------------------------------------
			// Ynew = conj(transpose(W_S)) * Yold * transpose(W_E)
			// ---------------------------------------------------
			//
			// Method 1:
			// (1) Ytemp = conj( transpose(W_S)) * Yold
			// (2) Ynew = Ytemp * transpose(W_E)
			//
			//  need   2.0 * nrow_W_S * ncol_W_S * ncol_Yold +
			//         2.0 * nrow_Ytemp * ncol_Ytemp *  nrow_W_E
			// Method 2:
			// (1) Ytemp = Yold * transpose(W_E)
			// (2) Ynew = conj(transpose(W_S)) * Ytemp
			//
			// need  2.0 * nrow_Yold * ncol_Yold  * nrow_W_E +
			//       2.0 * nrow_W_S * ncol_W_S * ncol_Ytemp
			// ---------------------------------------------------


			// ---------------------------------------
			// Method 1:
			// (1) Ytemp = conj( transpose(W_S)) * Yold
			// (2) Ynew = Ytemp * transpose(W_E)
			// ---------------------------------------
			nrow_Ytemp = ncol_W_S;
			ncol_Ytemp = ncol_Yold;
			const RealType flops_method_1 = 2.0 * nrow_W_S * ncol_W_S * ncol_Yold +
			        2.0 * nrow_Ytemp * ncol_Ytemp * nrow_W_E;

			// -----------------------
			// Method 2:
			// (1) Ytemp = Yold * transpose(W_E)
			// (2) Ynew = conj(transpose(W_S)) * Ytemp
			// -----------------------
			nrow_Ytemp = nrow_Yold;
			ncol_Ytemp = nrow_W_E;
			const RealType flops_method_2 = 2.0 * nrow_Yold * ncol_Yold * nrow_W_E +
			        2.0 * nrow_W_S * ncol_W_S * ncol_Ytemp;

			const bool use_method_1 = (flops_method_1 <= flops_method_2);

#if 0
			if (idebug >= 1) {
				std::cout << "WftAccelBlocks.h:360: "
				          << " use_method_1=" << use_method_1
				          << std::scientific
				          << " flops_method_1=" << flops_method_1
				          << " flops_method_2=" << flops_method_2
				          << std::defaultfloat
				          << "\n";

			}
#endif

			nrow_Ytemp = (use_method_1) ? ncol_W_S : nrow_Yold;
			ncol_Ytemp = (use_method_1) ? ncol_Yold : nrow_W_E;

			MatrixType tmp(nrow_Ytemp, ncol_Ytemp);
			tmp.setTo(0.0);

			ComplexOrRealType *Ytemp = &(tmp(0,0));
			const int ldYtemp = nrow_Ytemp;

			if (use_method_1) {
				// ---------------------------
				// Method 1:
				// (1) Ytemp = conj( transpose(W_S)) * Yold
				// (2) Ynew = Ytemp * transpose(W_E)
				// ---------------------------


				// -----------------------------------------
				// (1) Ytemp = conj( transpose(W_S)) * Yold
				// -----------------------------------------
				{
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					const int mm = nrow_Ytemp;
					const int nn = ncol_Ytemp;
					const int kk = nrow_Yold;

					psimag::BLAS::GEMM('C', 'N',
					                   mm,nn,kk,
					                   alpha, W_S, ldW_S, Yold, ldYold,
					                   beta,  Ytemp, ldYtemp );

				}

				// ----------------------------------
				// (2) Ynew = Ytemp * transpose(W_E)
				// ----------------------------------
				{
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					const int mm = nrow_Ynew;
					const int nn = ncol_Ynew;
					const int kk = ncol_Ytemp;

					psimag::BLAS::GEMM('N', 'T',
					                   mm,nn,kk,
					                   alpha, Ytemp, ldYtemp, W_E, ldW_E,
					                   beta, Ynew, ldYnew );

				}



			}
			else {
				// --------------------------------------
				// Method 2:
				// (1) Ytemp = Yold * transpose(W_E)
				// (2) Ynew = conj(transpose(W_S)) * Ytemp
				// --------------------------------------

				// -------------------------------
				// Ytemp = Yold * transpose( W_E )
				// -------------------------------

				{
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					const int mm = nrow_Ytemp;
					const int nn = ncol_Ytemp;
					const int kk = ncol_Yold;

					psimag::BLAS::GEMM('N', 'T',
					                   mm, nn, kk,
					                   alpha, Yold, ldYold, W_E,  ldW_E,
					                   beta, Ytemp, ldYtemp );
				}

				// ------------------------------------
				// Ynew = conj(transpose(W_S)) * Ytemp
				// ------------------------------------

				{
					const ComplexOrRealType alpha = d_one;
					const ComplexOrRealType beta = d_zero;

					const int mm = nrow_Ynew;
					const int nn = ncol_Ynew;
					const int kk = nrow_Ytemp;

					psimag::BLAS::GEMM('C', 'N',
					                   mm, nn, kk,
					                   alpha, W_S, ldW_S, Ytemp, ldYtemp,
					                   beta, Ynew, ldYnew);

				}
			};
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
		const ProgramGlobals::SysOrEnvEnum sysOrEnv_;
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
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toDense(ws);

		MatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toDense(we);

		SizeType i2psize = ws.cols();
		SizeType jp2size = we.rows();

		VectorMatrixType psi(volumeOfNk);
		for (SizeType kp = 0; kp < volumeOfNk; ++kp) {
			psi[kp].resize(i2psize, jp2size);
			psi[kp].setTo(0.0);
		}

		environPreparePsi(psi, psiSrc, i0src, volumeOfNk);

		VectorMatrixType result(volumeOfNk);

		SizeType threads = std::min(volumeOfNk,
		                            PsimagLite::Concurrency::codeSectionParams.npthreads);
		typedef PsimagLite::Parallelizer<ParallelWftInBlocks> ParallelizerType;
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerType threadedWft(codeSectionParams);

		ParallelWftInBlocks helperWft(result,
		                              psi,
		                              ws,
		                              we,
		                              volumeOfNk,
		                              ProgramGlobals::SysOrEnvEnum::ENVIRON,
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
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::SYSTEM).toDense(ws);

		MatrixType we;
		dmrgWaveStruct_.getTransform(ProgramGlobals::SysOrEnvEnum::ENVIRON).toDense(we);

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
		                              ProgramGlobals::SysOrEnvEnum::SYSTEM,
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
