#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include "Vector.h"
#include <numeric>
#include "BLAS.h"

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm2 {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef typename MatrixDenseOrSparseType::MatrixType MatrixType;
	typedef long int IntegerType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<char>::Type VectorCharType;
	typedef typename PsimagLite::Vector<ComplexOrRealType*>::Type VectorStarType;
	typedef typename PsimagLite::Vector<const ComplexOrRealType*>::Type VectorConstStarType;

	static const int ialign_ = 32;
	static const int idebug_ = 0; // set to 0 until it gives correct results

	class MatrixOrVector {

	public:

		MatrixOrVector(const VectorType& v) : v_(&v), m_(0) {}

		MatrixOrVector(const MatrixType& m) : v_(0), m_(&m) {}

		const ComplexOrRealType* operator()(SizeType first) const
		{
			if (v_) {
				assert(!m_);
				return &(v_->operator[](first));
			}

			assert(m_);
			return &(m_->operator()(first, 0));
		}

		const ComplexOrRealType& operator()(SizeType first,
		                                    SizeType row,
		                                    SizeType col,
		                                    SizeType ld) const
		{
			if (v_) {
				assert(!m_);
				SizeType offset = first + row + col*ld;
				assert(offset < v_->size());
				return v_->operator[](offset);
			}

			assert(m_);
			return m_->operator()(row + first, col);
		}

	private:

		const VectorType* v_;
		const MatrixType* m_;
	};

public:

	BatchedGemm2(const InitKronType& initKron) : initKron_(initKron)
	{
		if (!enabled()) return;

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
	}

	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		if (!enabled())
			err("BatchedGemm::matrixVector called but BatchedGemm not enabled\n");

		RealType gflops1 = 0.0;
		RealType gflops2 = 0.0;
		RealType time1stVbatch = 0.0;
		RealType time2ndVbatch = 0.0;
		/*
 ------------------
 compute  Y = H * X
 ------------------
*/
		int leftMaxStates  = initKron_.lrs(InitKronType::NEW).left().size();
		int rightMaxStates = initKron_.lrs(InitKronType::NEW).right().size();

		SizeType npatches = initKron_.numberOfPatches(InitKronType::OLD);
		SizeType noperator = initKron_.connections();
		VectorSizeType leftPatchSize(npatches);
		VectorSizeType rightPatchSize(npatches);

		for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType igroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::LEFT)[ipatch];

			int L1 = initKron_.lrs(InitKronType::NEW).left().partition(igroup);
			int L2 = initKron_.lrs(InitKronType::NEW).left().partition(igroup + 1);

			leftPatchSize[ipatch] =  L2 - L1;
		}

		for(SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			SizeType igroup = initKron_.patch(InitKronType::NEW,
			                                  GenIjPatchType::RIGHT)[ipatch];
			int  R1 = initKron_.lrs(InitKronType::NEW).right().partition(igroup);
			int  R2 = initKron_.lrs(InitKronType::NEW).right().partition(igroup + 1);
			rightPatchSize[ipatch] = R2 - R1;
		}

		int nrowA = leftMaxStates;
		int ncolA = nrowA;
		int nrowB = rightMaxStates;
		int ncolB = nrowB;
		int nrowBX = nrowB;
		int ncolBX = ncolA * noperator;
		int ldBX = ialign_ * iceil(nrowBX, ialign_);
		MatrixType BX(ldBX,  ncolA*noperator);

		MatrixOrVector mOvVin(vin);

		time1stVbatch = -dmrgGetWtime();

		for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
			long j1 = initKron_.offsetForPatches(InitKronType::NEW, jpatch);
			int nrowX = rightPatchSize[jpatch];
			assert(initKron_.offsetForPatches(InitKronType::NEW, jpatch + 1) - j1 ==
			       nrowX * leftPatchSize[jpatch]);

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
				int mm = nrowBX;
				int nn = L2 - L1;
				int kk = R2 - R1;

				gflops1 += ((2.0*mm)*nn)*kk;

				almostGemm('N',
				           mm,
				           nn,
				           kk,
				           Bbatch_,
				           offsetB + R1,
				           j1,
				           ldXJ,
				           mOvVin,
				           &(BX(0, offsetBX + L1)),
				           ldBX);

				gflops1 = gflops1/(1000.0*1000.0*1000.0);
			}
		}

		time1stVbatch += dmrgGetWtime();
		/*
 -------------------------------------------------
 perform computations with  Y += (BX)*transpose(A)
 -------------------------------------------------
*/
		time2ndVbatch = -dmrgGetWtime();

		MatrixOrVector mOvAbatch(Abatch_);

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

			assert(R2 - R1 == rightPatchSize[ipatch] &&
			       L2 - L1 == leftPatchSize[ipatch]);

			assert(static_cast<SizeType>(i1) < vout.size());
			ComplexOrRealType *YI = &(vout[i1]);
			int nrowYI = R2 - R1;
			int ldYI = nrowYI;
			int ncolYI = L2 - L1;
			assert(initKron_.offsetForPatches(InitKronType::NEW, ipatch + 1) - i1 ==
			       nrowYI * ncolYI);

			/*
		--------------------------------------------------------------------
		YI(1:(R2-R1+1),1:(L2-L1+1)) = BX( R1:R2,1:ncolBX) *
										 transpose( Abatch( L1:L2,1:ncolBX) );
		--------------------------------------------------------------------
	  */
			int mm = nrowYI;
			int nn = ncolYI;
			int kk = ncolBX;

			gflops2 += ((2.0*mm)*nn)*kk;

			almostGemm('T',
			           mm,
			           nn,
			           kk,
			           BX,
			           R1,
			           L1,
			           Abatch_.rows(),
			           mOvAbatch,
			           YI,
			           ldYI);
		}

		time2ndVbatch += dmrgGetWtime();
		gflops2 /= (1000.0*1000.0*1000.0);

		if (time1stVbatch < 1e-6 || time2ndVbatch < 1e-6)
			return;
		std::cerr<<"1st vbatch "<<gflops1/time1stVbatch;
		std::cerr<<" gflops/sec (gflops1="<<gflops1<<",time="<<time1stVbatch<<")\n";
		std::cerr<<"2nd vbatch "<<gflops2/time2ndVbatch<<" gflops/sec (gflops2=";
		std::cerr<<gflops2<<",time="<<time2ndVbatch<<")\n";

		std::cerr<<"overall "<<(gflops1+gflops2)/(time1stVbatch + time2ndVbatch);
		std::cerr<<" gflops/sec\n";
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


	static RealType dmrgGetWtime()
	{
		return clock()/CLOCKS_PER_SEC;
	}

	static void almostGemm(char transb, int m, int n, int ktotal,
	                       const MatrixType& A, SizeType rowOrCol,
	                       SizeType firstB, int ldb, const MatrixOrVector& B,
	                       ComplexOrRealType* C, int ldc)
	{
		const ComplexOrRealType* aptr = (transb == 'N') ? &(A(0, rowOrCol)) : &(A(rowOrCol, 0));
		psimag::BLAS::GEMM('N',
		                   transb,
		                   m,
		                   n,
		                   ktotal,
		                   1.0,
		                   aptr,
		                   A.rows(),
		                   B(firstB),
		                   ldb,
		                   0.0,
		                   C,
		                   ldc);
		//		for (int j = 0; j < n; ++j) {
		//			for (int i = 0; i < m; ++i) {
		//				C[i + j*ldc] = 0.0;
		//				for (int k = 0; k < ktotal; ++k) {
		//					ComplexOrRealType opAik =  (transb == 'N') ? A(i, k + rowOrCol) :
		//					                                             A(i + rowOrCol, k);
		//					ComplexOrRealType opBkj = (transb == 'N') ? B(firstB, k, j, ldb) :
		//					                                            B(firstB, j, k, ldb);
		//					C[i + j*ldc] += opAik*opBkj;
		//				}
		//			}
		//		}
	}

	static void almostGemmOld(char transb, int m, int n, int ktotal,
	                          const MatrixType& A, SizeType rowOrCol,
	                          SizeType firstB, int ldb, const MatrixOrVector& B,
	                          ComplexOrRealType* C, int ldc)
	{
		for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
				C[i + j*ldc] = 0.0;
				for (int k = 0; k < ktotal; ++k) {
					ComplexOrRealType opAik =  (transb == 'N') ? A(i, k + rowOrCol) :
					                                             A(i + rowOrCol, k);
					ComplexOrRealType opBkj = (transb == 'N') ? B(firstB, k, j, ldb) :
					                                            B(firstB, j, k, ldb);
					C[i + j*ldc] += opAik*opBkj;
				}
			}
		}
	}

	const InitKronType& initKron_;
	bool enabled_;
	MatrixType Abatch_;
	MatrixType Bbatch_;
};
}
#endif // BATCHEDGEMM_H
