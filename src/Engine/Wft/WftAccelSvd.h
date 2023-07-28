#ifndef WFT_ACCEL_SVD_H
#define WFT_ACCEL_SVD_H
#include "BLAS.h"
#include "ProgramGlobals.h"
#include "Vector.h"
#include <limits>

namespace Dmrg
{

template <typename WaveFunctionTransfBaseType>
class WftAccelSvd
{

	typedef typename WaveFunctionTransfBaseType::DmrgWaveStructType DmrgWaveStructType;
	typedef typename DmrgWaveStructType::WaveStructSvdType WaveStructSvdType;
	typedef typename WaveFunctionTransfBaseType::WftOptionsType WftOptionsType;
	typedef typename WaveFunctionTransfBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename WaveFunctionTransfBaseType::VectorSizeType VectorSizeType;
	using OneSiteSpacesType = typename WaveFunctionTransfBaseType::OneSiteSpacesType;
	typedef typename DmrgWaveStructType::LeftRightSuperType LeftRightSuperType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef typename DmrgWaveStructType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::VectorQnType VectorQnType;
	typedef typename BasisWithOperatorsType::QnType QnType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename WaveFunctionTransfBaseType::PackIndicesType PackIndicesType;
	typedef typename DmrgWaveStructType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BlockDiagonalMatrixType::BuildingBlockType MatrixType;
	typedef typename PsimagLite::Vector<MatrixType>::Type VectorMatrixType;
	typedef typename PsimagLite::Vector<VectorType>::Type VectorVectorType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	class LoopOne
	{

	public:

		LoopOne(const MatrixType& u,
		    const MatrixType& vPrime,
		    const VectorQnType& q,
		    const MatrixType& uTildePrime,
		    const MatrixType& vTildePrimePrime,
		    const VectorQnType& qTilde)
		    : u_(u)
		    , vPrime_(vPrime)
		    , q_(q)
		    , uTildePrime_(uTildePrime)
		    , vTildePrimePrime_(vTildePrimePrime)
		    , qTilde_(qTilde)
		    , npatches_(std::max(u_.size(), uTildePrime_.size()))
		    , uFinal_(std::min(u_.size(), uTildePrime_.size()))
		{
			assert(u.size() == vPrime.size());
			assert(u.size() == q.size());
			assert(uTildePrime.size() == vTildePrimePrime.size());
			asserT(uTildePrime.size() == qTilde.size());
		}

		SizeType tasks() const { return npatches_; }

		void doTask(SizeType patchBig, SizeType)
		{
			const VectorQnType& qSmall = (q_.size() > qTilde_.size()) ? qTilde_ : q_;
			const VectorQnType& qBig = (q_.size() > qTilde_.size()) ? q_ : qTilde_;
			const QnType& q1 = qBig[patchBig];
			int patchSmall = indexOrMinusOne(qSmall, q1);
			if (patchSmall < 0)
				return;

			SizeType patch = (q_.size() > qTilde_.size()) ? patchBig : patchSmall;
			SizeType patchTilde = (q_.size() > qTilde_.size()) ? patchSmall : patchBig;

			uFinal_[patchSmall] = uTildePrime_[patchTilde] * u_[patch];
			vPrimeFinal_[patchSmall] = vPrime_[patchTilde] * vTildePrimePrime_[patch];
		}

		const VectorMatrixType& uFinal() const { return uFinal_; }

		const VectorMatrixType& vPrimeFinal() const { return vPrimeFinal_; }

		const VectorQnType& qns() const { return (q_.size() > qTilde_.size()) ? qTilde_ : q_; }

	private:

		const VectorMatrixType& u_;
		const VectorMatrixType& vPrime_;
		const VectorQnType& q_;
		const VectorMatrixType& uTildePrime_;
		const VectorMatrixType& vTildePrimePrime_;
		const VectorQnType& qTilde_;
		SizeType npatches_;
		VectorMatrixType uFinal_;
		VectorMatrixType vPrimeFinal_;
	}; // class LoopOne

	class LoopTwo
	{

	public:

		LoopTwo(const MatrixType& ufinal,
		    const MatrixType& vPrimeFinal,
		    const VectorQnType& qnsFinal,
		    const VectorType& d,
		    const VectorQnType& qnsOfD)
		    : uFinal_(ufinal)
		    , vPrimeFinal_(vPrimeFinal)
		    , qnsFinal_(qnsFinal)
		    , d_(d)
		    , qnsOfD_(qnsOfD)
		    , patches_(std::max(qnsFinal.size(), qnsOfD.size()))
		    , result_(std::min(qnsFinal.size(), qnsOfD.size()))
		{
			assert(ufinal.size() == vPrimeFinal.size());
			assert(ufinal.size() == qnsFinal.size());
			assert(d.size() == qnsOfD.size());
		}

		SizeType tasks() const { return patches_; }

		// U[patch] D V'[patch] *D[patch]
		void doTask(SizeType patchBig, SizeType)
		{
			const VectorQnType& qSmall = (qnsFinal_.size() > qnsOfD_.size()) ? qnsOfD_ : qnsFinal_;
			const VectorQnType& qBig = (qnsFinal_.size() > qnsOfD_.size()) ? qnsFinal_ : qnsOfD_;
			const QnType& q1 = qBig[patchBig];
			int patchSmall = indexOrMinusOne(qSmall, q1);
			if (patchSmall < 0)
				return;

			SizeType patch = (qnsFinal_.size() > qnsOfD_.size()) ? patchBig : patchSmall;
			SizeType patchOfD = (qnsFinal_.size() > qnsOfD_.size()) ? patchSmall : patchBig;

			const MatrixType& u = uFinal_[patch];
			const MatrixType& vprime = vPrimeFinal_[patch];
			const VectorType& d = d_[patchOfD];
			SizeType itotal = u.rows();
			SizeType ktotal = u.cols();
			assert(ktotal = d.size());
			SizeType jtotal = vprime.cols();
			assert(vprime.rows() == ktotal);

			for (SizeType i = 0; i < itotal; ++i) {
				for (SizeType k = 0; k < ktotal; ++k) {
					ComplexOrRealType prod = u(i, k) * d[k];
					for (SizeType j = 0; j < jtotal; ++j) {
						result_[patchSmall](i, j) += prod * vprime(k, j);
					}
				}
			}
		}

	private:

		const VectorMatrixType& uFinal_;
		const VectorMatrixType& vPrimeFinal_;
		const VectorQnType& qnsFinal_;
		const VectorVectorType& d_;
		const VectorQnType& qnsOfD_;
		SizeType patches_;
		VectorMatrixType result_;
	}; // class LoopTwo

public:

	WftAccelSvd(const DmrgWaveStructType& dmrgWaveStruct,
	    const WftOptionsType& wftOptions)
	    : dmrgWaveStruct_(dmrgWaveStruct)
	    , wftOptions_(wftOptions)
	{
	}

	void operator()(VectorWithOffsetType& psiDest,
	    SizeType iNew,
	    const VectorWithOffsetType& psiSrc,
	    SizeType iOld,
	    const LeftRightSuperType& lrs,
	    const OneSiteSpacesType& nk) const
	{
		//		typename ProgramGlobals::SysOrEnvEnum prevPart = (dir == ProgramGlobals::EXPAND_SYSTEM) ?
		//		            ProgramGlobals::SYSTEM : ProgramGlobals::ENVIRON;
		//		const WaveStructSvdType& wavePrev = dmrgWaveStruct_.getWave(prevPart);

		//		typename ProgramGlobals::SysOrEnvEnum oppoPart = (dir == ProgramGlobals::EXPAND_SYSTEM) ?
		//		            ProgramGlobals::ENVIRON : ProgramGlobals::SYSTEM;
		//		const WaveStructSvdType& waveOld = dmrgWaveStruct_.getWave(oppoPart);

		//		internal(waveOld.u(),
		//		         waveOld.vts(),
		//		         waveOld.qns(),
		//		         wavePrev.u(),
		//		         wavePrev.vts(),
		//		         wavePrev.qns(),
		//		         wavePrev.s());
	}

private:

	void internal(const VectorMatrixType& uVeryOld,
	    const VectorMatrixType& vPrimeVeryOld,
	    const VectorQnType& qnsVeryOld,
	    const VectorMatrixType& uPrevious,
	    const VectorMatrixType& vPrimePrevious,
	    const VectorQnType& qnsPrevious,
	    const VectorVectorType& sPrevious)
	{

		err("WftAccelSvd: Not ready yet\n");
		VectorMatrixType uPreviousPinv;
		pinv(uPreviousPinv, uPrevious);

		VectorMatrixType vPrimePreviousPinv;
		pinv(vPrimePreviousPinv, vPrimePrevious);
		LoopOne loopOne(uVeryOld,
		    vPrimeVeryOld,
		    qnsVeryOld,
		    uPreviousPinv,
		    vPrimePreviousPinv,
		    qnsPrevious);
		typedef PsimagLite::Parallelizer<LoopOne> ParallelizerOneType;
		SizeType threads = std::min(std::max(qnsVeryOld.size(), qnsPrevious.size()),
		    PsimagLite::Concurrency::codeSectionParams.npthreads);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		ParallelizerOneType threadOne(codeSectionParams);
		threadOne.loopCreate(loopOne);

		typedef PsimagLite::Parallelizer<LoopTwo> ParallelizerTwoType;
		LoopTwo loopTwo(loopOne.uFinal(),
		    loopOne.vPrimeFinal(),
		    loopOne.qns(),
		    sPrevious,
		    qnsPrevious);
		ParallelizerTwoType threadTwo(codeSectionParams);
		threadTwo.loopCreate(loopTwo);
	}

	static void pinv(VectorMatrixType& dest, const VectorMatrixType& src)
	{
		SizeType n = src.size();
		dest.resize(n);
		for (SizeType i = 0; i < n; ++i)
			pinvOne(dest[i], src[i]);
	}

	static void pinvOne(MatrixType& dest, const MatrixType& src)
	{
		dest = src;
		VectorRealType s;
		MatrixType vt;
		svd('A', dest, s, vt);
		SizeType rows = dest.rows();
		SizeType cols = vt.cols();
		RealType epsilon = *std::max_element(s.begin(), s.end()) * std::max(rows, cols) * std::numeric_limits<RealType>::epsilon();

		SizeType n = biggerThanEpsilon(s, epsilon);

		MatrixType uoneTranspose(n, rows);
		MatrixType voneTimeStoTheMinusOne(cols, n);
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType r = 0; r < rows; ++r)
				uoneTranspose(i, r) = dest(r, i);
			for (SizeType c = 0; c < cols; ++c)
				voneTimeStoTheMinusOne(c, i) = vt(i, c) / s[i];
		}

		dest.resize(cols, rows);
		ComplexOrRealType alpha = 1;
		ComplexOrRealType beta = 0;
		psimag::BLAS::GEMM('N',
		    'N',
		    cols,
		    rows,
		    n,
		    alpha,
		    &(voneTimeStoTheMinusOne(0, 0)),
		    cols,
		    &(uoneTranspose(0, 0)),
		    n,
		    beta,
		    dest,
		    cols);
	}

	static SizeType biggerThanEpsilon(const VectorRealType& v, RealType epsilon)
	{
		SizeType n = v.size();
		for (SizeType i = 0; i < n; ++i)
			if (std::abs(v[i]) < epsilon)
				return i;
		return n;
	}

	const DmrgWaveStructType& dmrgWaveStruct_;
	const WftOptionsType& wftOptions_;
}; // class WftAccelSvd
}
#endif // WFT_ACCEL_SVD_H
