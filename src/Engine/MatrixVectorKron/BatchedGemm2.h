#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include "Vector.h"
#include <numeric>
#include "BLAS.h"
#include "ProgressIndicator.h"

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm2 {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef long int IntegerType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<char>::Type VectorCharType;
	typedef typename PsimagLite::Vector<ComplexOrRealType*>::Type VectorStarType;
	typedef typename PsimagLite::Vector<const ComplexOrRealType*>::Type VectorConstStarType;

	static const int ialign_ = 32;
	static const int idebug_ = 0; // set to 0 until it gives correct results

public:

	BatchedGemm2(const InitKronType& initKron)
	    : initKron_(initKron), progress_("BatchedGemm")
	{
		if (!enabled()) return;

		{
			PsimagLite::OstringStream msg;
			msg<<"Constructing...";
			progress_.printline(msg,std::cout);
		}

		SizeType npatches = initKron_.numberOfPatches(InitKronType::OLD);
		SizeType noperator = initKron_.connections();

		SizeType leftMaxState = initKron_.lrs(InitKronType::NEW).left().size();
		SizeType rightMaxState = initKron_.lrs(InitKronType::NEW).right().size();

		int nrowAbatch = leftMaxState;
		int ncolAbatch = leftMaxState * noperator;

		int nrowBbatch = rightMaxState;
		int ncolBbatch = rightMaxState * noperator;

		int ldAbatch = ialign_ * iceil(nrowAbatch, ialign_ );
		int ldBbatch = ialign_ * iceil(nrowBbatch, ialign_ );

		Abatch_.resize(ldAbatch, ncolAbatch);
		Bbatch_.resize(ldBbatch, ncolBbatch);

		/*
  -------------------------
  fill in Abatch and Bbatch
  -------------------------
  */

		assert(ldAbatch * leftMaxState * noperator >= 1);
		assert(ldBbatch * rightMaxState * noperator >= 1);

		for (SizeType ioperator = 0; ioperator < noperator; ++ioperator) {
			const ArrayOfMatStructType& xiStruct = initKron_.xc(ioperator);
			for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
				for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {

					const MatrixType& Asrc =  xiStruct(ipatch,jpatch).dense();
					SizeType igroup = initKron_.patch(InitKronType::NEW,
					                                  GenIjPatchType::LEFT)[ipatch];
					SizeType jgroup = initKron_.patch(InitKronType::NEW,
					                                  GenIjPatchType::LEFT)[jpatch];
					int ia = initKron_.lrs(InitKronType::NEW).left().partition(igroup);
					int ja = initKron_.lrs(InitKronType::NEW).left().partition(jgroup);

					mylacpy(Asrc, Abatch_, ia, ja + ioperator*leftMaxState);
				}
			}
		}

		for (SizeType ioperator = 0; ioperator < noperator; ++ioperator) {
			const ArrayOfMatStructType& yiStruct = initKron_.yc(ioperator);
			for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
				for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {

					const MatrixType& Bsrc =  yiStruct(ipatch,jpatch).dense();
					SizeType igroup = initKron_.patch(InitKronType::NEW,
					                                  GenIjPatchType::RIGHT)[ipatch];
					SizeType jgroup = initKron_.patch(InitKronType::NEW,
					                                  GenIjPatchType::RIGHT)[jpatch];
					int ib = initKron_.lrs(InitKronType::NEW).right().partition(igroup);
					int jb = initKron_.lrs(InitKronType::NEW).right().partition(jgroup);

					mylacpy(Bsrc, Bbatch_, ib, jb + ioperator*rightMaxState);
				}
			}
		}

		leftPatchSize_.resize(npatches, 0);
		rightPatchSize_.resize(npatches, 0);

		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType igroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::LEFT)[ipatch];

			int L1 = initKron_.lrs(InitKronType::NEW).left().partition(igroup);
			int L2 = initKron_.lrs(InitKronType::NEW).left().partition(igroup + 1);

			leftPatchSize_[ipatch] =  L2 - L1;
		}

		for(SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType igroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::RIGHT)[ipatch];
			int  R1 = initKron_.lrs(InitKronType::NEW).right().partition(igroup);
			int  R2 = initKron_.lrs(InitKronType::NEW).right().partition(igroup + 1);
			rightPatchSize_[ipatch] = R2 - R1;
		}

		int leftMaxStates  = initKron_.lrs(InitKronType::NEW).left().size();
		int rightMaxStates = initKron_.lrs(InitKronType::NEW).right().size();
		int nrowA = leftMaxStates;
		int ncolA = nrowA;
		int nrowB = rightMaxStates;
		int nrowBX = nrowB;
		int ldBX = ialign_ * iceil(nrowBX, ialign_);
		BX_.resize(ldBX,  ncolA*noperator);

		{
			PsimagLite::OstringStream msg;
			msg<<"Construction done.";
			progress_.printline(msg,std::cout);
		}
	}

	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		if (!enabled())
			err("BatchedGemm::matrixVector called but BatchedGemm not enabled\n");

		/*
 ------------------
 compute  Y = H * X
 ------------------
*/
		int leftMaxStates  = initKron_.lrs(InitKronType::NEW).left().size();
		int rightMaxStates = initKron_.lrs(InitKronType::NEW).right().size();
		SizeType npatches = initKron_.numberOfPatches(InitKronType::OLD);
		SizeType noperator = initKron_.connections();
		int nrowA = leftMaxStates;
		int ncolA = nrowA;
		int nrowB = rightMaxStates;
		int ncolB = nrowB;
		int nrowBX = nrowB;
		int ncolBX = ncolA * noperator;
		int ldBX = ialign_ * iceil(nrowBX, ialign_);
		BX_.setTo(0.0);

		for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
			long j1 = initKron_.offsetForPatches(InitKronType::NEW, jpatch);
			int nrowX = rightPatchSize_[jpatch];
			assert(initKron_.offsetForPatches(InitKronType::NEW, jpatch + 1) - j1 ==
			       nrowX * leftPatchSize_[jpatch]);

			/*
	 --------------------------------------
	 XJ = reshape( X(j1:j2), nrowX, ncolX )
	 --------------------------------------
	 */
			assert(static_cast<SizeType>(j1) < vin.size());
			int ldXJ = nrowX;

			SizeType jgroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::RIGHT)[jpatch];
			int R1 = initKron_.lrs(InitKronType::NEW).right().partition(jgroup);
			int R2 = initKron_.lrs(InitKronType::NEW).right().partition(jgroup + 1);

			SizeType igroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::LEFT)[jpatch];
			int L1 = initKron_.lrs(InitKronType::NEW).left().partition(igroup);
			int L2 = initKron_.lrs(InitKronType::NEW).left().partition(igroup + 1);

			assert(static_cast<SizeType>(j1 + R2 - R1 - 1 + (L2 - L1 - 1)*nrowX) <
			       vin.size());
			/*
	 -------------------------------
	 independent DGEMM in same group
	 -------------------------------
	 */
			for (SizeType k = 0; k < noperator; ++k) {
				int offsetB = k*ncolB;
				int offsetBX = k*ncolA;
				/*
		------------------------------------------------------------------------
		BX(1:nrowBX, offsetBX + (L1:L2)) = Bbatch(1:nrowBX, offsetB + (R1:R2) ) *
											 XJ( 1:(R2-R1+1), 1:(L2-L1+1));
		------------------------------------------------------------------------
		*/
				psimag::BLAS::GEMM('N',
				                   'N',
				                   nrowBX,
				                   L2 - L1,
				                   R2 - R1,
				                   1.0,
				                   &(Bbatch_(0, offsetB + R1)),
				                   Bbatch_.rows(),
				                   &(vin[j1]),
				                   ldXJ,
				                   0.0,
				                   &(BX_(0, offsetBX + L1)),
				                   ldBX);
			}
		}

		/*
 -------------------------------------------------
 perform computations with  Y += (BX)*transpose(A)
 -------------------------------------------------
*/
		for(SizeType ipatch = 0; ipatch < npatches; ++ipatch) {

			long i1 = initKron_.offsetForPatches(InitKronType::NEW, ipatch);


			SizeType jgroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::RIGHT)[ipatch];
			SizeType R1 = initKron_.lrs(InitKronType::NEW).right().partition(jgroup);
			SizeType R2 = initKron_.lrs(InitKronType::NEW).right().partition(jgroup + 1);

			SizeType igroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::LEFT)[ipatch];
			SizeType L1 = initKron_.lrs(InitKronType::NEW).left().partition(igroup);
			SizeType L2 = initKron_.lrs(InitKronType::NEW).left().partition(igroup + 1);

			assert(R2 - R1 == rightPatchSize_[ipatch] &&
			       L2 - L1 == leftPatchSize_[ipatch]);

			assert(static_cast<SizeType>(i1) < vout.size());
			ComplexOrRealType *YI = &(vout[i1]);
			int nrowYI = R2 - R1;
			int ldYI = nrowYI;
			int ncolYI = L2 - L1;
			assert(static_cast<int>(initKron_.offsetForPatches(InitKronType::NEW, ipatch + 1) - i1) ==
			       nrowYI * ncolYI);

			/*
		--------------------------------------------------------------------
		YI(1:(R2-R1+1),1:(L2-L1+1)) = BX( R1:R2,1:ncolBX) *
										 transpose( Abatch( L1:L2,1:ncolBX) );
		--------------------------------------------------------------------
	  */
			psimag::BLAS::GEMM('N',
			                   'T',
			                   nrowYI,
			                   ncolYI,
			                   ncolBX,
			                   1.0,
			                   &(BX_(R1, 0)),
			                   BX_.rows(),
			                   &(Abatch_(L1, 0)),
			                   Abatch_.rows(),
			                   0.0,
			                   YI,
			                   ldYI);
		}
	}

private:

	static int iceil(int x, int n)
	{
		return (x + n - 1)/n;
	}

	static void mylacpy(const MatrixType& a,
	                    MatrixType& b,
	                    SizeType xstart,
	                    SizeType ystart)
	{
		int m = a.rows();
		int n = a.cols();
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < m; ++i)
				b(i + xstart, j + ystart) = a(i, j);
	}

	const InitKronType& initKron_;
	PsimagLite::ProgressIndicator progress_;
	MatrixType Abatch_;
	MatrixType Bbatch_;
	mutable MatrixType BX_;
	VectorSizeType leftPatchSize_;
	VectorSizeType rightPatchSize_;
};
}
#endif // BATCHEDGEMM_H
