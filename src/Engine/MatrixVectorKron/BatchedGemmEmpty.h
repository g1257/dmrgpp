#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include "Vector.h"

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm2 {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename VectorType::value_type ComplexOrRealType;
	typedef typename MatrixDenseOrSparseType::MatrixType MatrixType;
	typedef long int IntegerType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef PsimagLite::Vector<char>::Type VectorCharType;

public:

	BatchedGemm2(const InitKronType& initKron)
	    : initKron_(initKron),
	      xyPatchStart_(initKron_.numberOfPatches(InitKronType::OLD), 0)
	{
		if (!enabled()) return;

		SizeType npatches = initKron_.numberOfPatches(InitKronType::OLD);
		SizeType noperator = initKron_.connections();

		SizeType leftMaxState = 0;
		SizeType rightMaxState = 0;
		for(SizeType ipatch = 0; ipatch <  npatches; ++ipatch) {
			leftMaxState += initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
			rightMaxState += initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
		}

		int nrowAbatch = leftMaxState;
		int ncolAbatch = leftMaxState * noperator;

		int nrowBbatch = rightMaxState;
		int ncolBbatch = rightMaxState * noperator;

		const int ialign = 32;
		int ldAbatch = ialign * iceil(nrowAbatch, ialign );
		int ldBbatch = ialign * iceil(nrowBbatch, ialign );

		MatrixType Abatch(ldAbatch, ncolAbatch);
		MatrixType Bbatch(ldBbatch, ncolBbatch);

		for (SizeType ipatch = 0; ipatch < npatches - 1; ++ipatch) {
			int nrowX = initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
			int ncolX = initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
			xyPatchStart_[ipatch + 1] = xyPatchStart_[ipatch] + nrowX*ncolX;
		}

		/*
  -------------------------
  fill in Abatch and Bbatch
  -------------------------
  */

		const size_t sizeAbatch = ldAbatch * leftMaxState * noperator;
		const size_t sizeBbatch = ldBbatch * rightMaxState * noperator;

		assert(sizeAbatch >= 1);
		assert(sizeBbatch >= 1);

		//#ifdef USE_GETSET
		//		ComplexOrRealType *hAbatch_ = (ComplexOrRealType *) malloc( sizeof(ComplexOrRealType) * sizeAbatch );
		//#else
		//		ComplexOrRealType *hAbatch_ = Abatch_;
		//#endif

		for (SizeType ioperator = 0; ioperator < noperator; ++ioperator) {
			const ArrayOfMatStructType& xiStruct = initKron_.xc(ioperator);
			for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
				for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {

					const MatrixType& Asrc =  xiStruct(ipatch,jpatch).dense();

					int nrowAsrc = initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1)
					        - initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
					int ncolAsrc = initKron_.lrs(InitKronType::NEW).left().partition(jpatch + 1)
					        - initKron_.lrs(InitKronType::NEW).left().partition(jpatch);

					int ia = initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
					int ja = initKron_.lrs(InitKronType::NEW).left().partition(jpatch);

					ComplexOrRealType *Adest = &(Abatch(ia, ja + ioperator*leftMaxState));

					SizeType tmp = nrowAsrc - 1 + (ncolAsrc - 1)*ldAbatch;
					assert(tmp < Abatch.rows()*Abatch.cols());
					std::cerr<<"HERE "<<ioperator<<" "<<jpatch<<" "<<ipatch<<"\n";
					mylacpy(nrowAsrc, ncolAsrc, Asrc, Adest, ldAbatch);
				}
			}
		}

		// USE_GETSET block omitted (2 blocks)

		for (SizeType ioperator = 0; ioperator < noperator; ++ioperator) {
			const ArrayOfMatStructType& yiStruct = initKron_.yc(ioperator);
			for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
				for (SizeType ipatch = 0; ipatch < npatches; ++ipatch) {

					const MatrixType& Bsrc =  yiStruct(ipatch,jpatch).dense();

					int nrowBsrc = initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1)
					        - initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
					int ncolBsrc = initKron_.lrs(InitKronType::NEW).right().partition(jpatch + 1)
					        - initKron_.lrs(InitKronType::NEW).right().partition(jpatch);

					int ib = initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
					int jb = initKron_.lrs(InitKronType::NEW).right().partition(jpatch);

					ComplexOrRealType* Bdest = &(Bbatch(ib,jb + ioperator*rightMaxState));
					SizeType tmp = nrowBsrc - 1 + (ncolBsrc - 1)*ldBbatch;
						assert(tmp < Bbatch.rows()*Bbatch.cols());
					mylacpy(nrowBsrc, ncolBsrc, Bsrc, Bdest, ldBbatch);
				}
			}
		}

		// USE_GETSET block here omitted

//		*pAbatch = Abatch_;
//		*pBbatch = Bbatch_;
//		*pleftPatchStart = leftPatchStart_;
//		*prightPatchStart = rightPatchStart_;
//		*pxyPatchStart = xyPatchStart_;
	}


	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType&, const VectorType&) const
	{
		if (!enabled())
			err("BatchedGemm::matrixVector called but BatchedGemm not enabled\n");

		const int ialign = 32;
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
			int L1 = initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
			int L2 = initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1);

			leftPatchSize[ipatch] =  L2 - L1;
		}

		for(SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			int  R1 = initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
			int  R2 = initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1);
			rightPatchSize[ipatch] = R2 - R1;
		}

		int ngroups = npatches;
		int ngroupsDim = ialign * iceil(ngroups, ialign);
		int batchSize = ngroups * noperator;
		int batchSizeDim = ialign * iceil(batchSize, ialign);

		VectorType alphaArray(ngroupsDim);
		VectorType betaArray(ngroupsDim);
//		ComplexOrRealType *aArray_[batchSizeDim];
//		ComplexOrRealType *bArray_[batchSizeDim];
//		ComplexOrRealType *cArray_[batchSizeDim];

		VectorSizeType mArray(ngroupsDim);
		VectorSizeType nArray(ngroupsDim);
		VectorSizeType kArray(ngroupsDim);
		VectorSizeType groupSize(ngroupsDim);
		VectorSizeType ldaArray(batchSizeDim);
		VectorSizeType ldbArray(batchSizeDim);
		VectorSizeType ldcArray(batchSizeDim);

		VectorCharType transaArray(ngroupsDim);
		VectorCharType transbArray(ngroupsDim);

		int nrowA = leftMaxStates;
		int ncolA = nrowA;
		int nrowB = rightMaxStates;
		int ncolB = nrowB;

		int nrowBX = nrowB;
		int ncolBX = ncolA * noperator;
		int ldBX = ialign * iceil(nrowBX, ialign);

// ComplexOrRealType *BX_ = (ComplexOrRealType *) dmrg_malloc( (sizeof(ComplexOrRealType) * ldBX) * (ncolA * noperator) );
// assert( BX_ != NULL );

//#define BX(i,j) BX_[ indx2f(i,j,ldBX) ]
		SizeType idx = 0;
		for (SizeType jpatch = 0; jpatch < npatches; ++jpatch) {
			int igroup = jpatch;
			long j1 = xyPatchStart_[jpatch];
			long j2 = xyPatchStart_[jpatch + 1];
			int nrowX = rightPatchSize[jpatch];
			int ncolX = leftPatchSize[jpatch];
			assert(j2 - j1 == nrowX * ncolX);

    /*
     --------------------------------------
     XJ = reshape( X(j1:j2), nrowX, ncolX )
     --------------------------------------
     */
			ComplexOrRealType *XJ = &( X(j1) );
			int ldXJ = nrowX;

			int R1 = initKron_.lrs(InitKronType::NEW).right().partition(jpatch);
			int R2 = initKron_.lrs(InitKronType::NEW).right().partition(jpatch + 1);
			int L1 = initKron_.lrs(InitKronType::NEW).left().partition(jpatch);
			int L2 = initKron_.lrs(InitKronType::NEW).left().partition(jpatch + 1);
			int kmax = noperator;

    /*
     -------------------------------
     independent DGEMM in same group
     -------------------------------
     */
			assert(igroup < groupSize.size());
			groupSize[igroup] = kmax;
			for (SizeType k = 0; k < kmax; ++k) {
				int offsetB = (k - 1)*ncolB;
				int offsetBX = (k - 1)*ncolA;

        /*
        ------------------------------------------------------------------------
        BX(1:nrowBX, offsetBX + (L1:L2)) = Bbatch(1:nrowBX, offsetB + (R1:R2) ) *
                                             XJ( 1:(R2-R1+1), 1:(L2-L1+1));
        ------------------------------------------------------------------------
        */
				transaArray[igroup] = 'N';
				transbArray[igroup] = 'N';
				int mm = nrowBX;
				int nn = L2 - L1;
				int kk = R2 - R1;
				mArray[igroup] = mm;
				nArray[igroup] = nn;
				kArray[igroup] = kk;

				gflops1 += ((2.0*mm)*nn)*kk;

				alphaArray[igroup] = 1.0;
				betaArray[igroup] = 0.0;

				cArray(idx) = &(BX(1,offsetBX+L1));
				ldcArray[igroup] = ldBX;

				aArray(idx) = &(Bbatch(1,offsetB+R1));
				ldaArray[igroup] = ld_Bbatch;

				bArray(idx) = XJ;
				ldbArray[igroup] = ldXJ;
				++idx;
			}
		}
   /*
    ------------------
    first vbatch DGEMM
    ------------------
    */
		time1stVbatch = -dmrg_get_wtime();
		dmrg_Xgemm_vbatch( transaArray, transbArray,
		                   mArray, nArray, kArray,
		                   alphaArray_,  aArray_, ldaArray, bArray_, ldbArray,
		                   betaArray,   cArray_, ldcArray,
		                   ngroups, groupSize );
		time1stVbatch += dmrg_get_wtime();
		gflops1 = gflops1/(1000.0*1000.0*1000.0);

/*
 -------------------------------------------------
 perform computations with  Y += (BX)*transpose(A)
 -------------------------------------------------
*/
		for(SizeType ipatch = 0; ipatch < npatches; ++ipatch) {
			int igroup = ipatch;

			long i1 = xyPatchStart_[ipatch];
			long i2 = xyPatchStart_[ipatch + 1];

			int R1 = initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
			int R2 = initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1);

			int L1 = initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
			int L2 = initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1);

			assert(R2 - R1 == rightPatchSize[ipatch] &&
			       L2 - L1 == leftPatchSize[ipatch]);

			//     ComplexOrRealType *YI = &(Y(i1));
			int nrowYI = R2 - R1;
			int ldYI = nrowYI;
			int ncolYI = L2 - L1;
			assert(i2 - i1 == nrowYI * ncolYI);

     /*
        --------------------------------------------------------------------
        YI(1:(R2-R1+1),1:(L2-L1+1)) = BX( R1:R2,1:ncolBX) *
                                         transpose( Abatch( L1:L2,1:ncolBX) );
        --------------------------------------------------------------------
      */
			groupSize[igroup] = 1;
			transaArray[igroup] = 'N';
			transbArray[igroup] = 'T';
			int mm = nrowYI;
			int nn = ncolYI;
			int kk = ncolBX;
			mArray[igroup] = mm;
			nArray[igroup] = nn;
			kArray[igroup] = kk;
			gflops2 += ((2.0*mm)*nn)*kk;
			alphaArray[igroup] = 1.0;
			betaArray[igroup] = 0.0;
			aArray[igroup] =  &(BX(R1,1));
			ldaArray[igroup] = ldBX;
			bArray[igroup] = &(Abatch(L1,1));
			ldbArray[igroup] = ld_Abatch;
			cArray[igroup] = YI;
			ldcArray[igroup] = ldYI;
		}

		ngroups = npatches;

   /*
    ------------------
    second vbatch DGEMM
    ------------------
	*/
		time2ndVbatch = -dmrg_get_wtime();
		dmrg_Xgemm_vbatch( transaArray, transbArray,
		                   mArray, nArray, kArray,
		                   alphaArray_,  aArray_, ldaArray, bArray_, ldbArray,
		                   betaArray,   cArray_, ldcArray,
		                   ngroups, groupSize );
		time2ndVbatch += dmrg_get_wtime();
		gflops2 = gflops2/(1000.0*1000.0*1000.0);

		std::cerr<<"1st vbatch "<<gflops1/time1stVbatch;
		std::cerr<<" gflops/sec (gflops1="<<gflops1<<",time="<<time1stVbatch<<")\n";
		std::cerr<<"2nd vbatch "<<gflops2/time2ndVbatch<<" gflops/sec (gflops2=";
		std::cerr<<gflops2<<",time="<<time2ndVbatch<<")\n";

		std::cerr<<"overall "<<(gflops1+gflops2)/(time1stVbatch + time2ndVbatch);
		std::cerr<<" gflops/sec\n";
// dmrg_free( (void *) BX_ );
	}

private:

	static int iceil(int x, int n)
	{
		return (x + n - 1)/n;
	}

	static void mylacpy(int m, int n, const MatrixType& a, ComplexOrRealType *b, int ldb)
	{
		for (int j = 0; j < n; ++j)
			for (int i = 0; i < m; ++i)
				b[i + j*(ldb)] = a(i, j);
	}

	const InitKronType& initKron_;
	bool enabled_;
	VectorSizeType xyPatchStart_;
};
}
#endif // BATCHEDGEMM_H
