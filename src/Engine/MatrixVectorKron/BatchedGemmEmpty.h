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
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	BatchedGemm2(const InitKronType& initKron) : initKron_(initKron)
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

		VectorSizeType leftPatchStart(npatches, 0);
		for (SizeType ipatch = 0; ipatch < npatches - 1; ++ipatch)
			leftPatchStart[ipatch + 1] = leftPatchStart[ipatch] +
			        initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).left().partition(ipatch);

		VectorSizeType rightPatchStart(npatches, 0);
		for (SizeType ipatch = 0; ipatch < npatches - 1; ++ipatch)
			rightPatchStart[ipatch + 1] = rightPatchStart[ipatch] +
			        initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).right().partition(ipatch);

		VectorSizeType xyPatchStart(npatches, 0);
		for (SizeType ipatch = 0; ipatch < npatches - 1; ++ipatch) {
			int nrowX = initKron_.lrs(InitKronType::NEW).right().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).right().partition(ipatch);
			int ncolX = initKron_.lrs(InitKronType::NEW).left().partition(ipatch + 1)
			        - initKron_.lrs(InitKronType::NEW).left().partition(ipatch);
			xyPatchStart[ipatch + 1] = xyPatchStart[ipatch] + nrowX*ncolX;
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
		//		FpType *hAbatch_ = (FpType *) malloc( sizeof(FpType) * sizeAbatch );
		//#else
		//		FpType *hAbatch_ = Abatch_;
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

					int ia = leftPatchStart[ipatch];
					int ja = leftPatchStart[jpatch];

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

					int ib = rightPatchStart[ipatch];
					int jb = rightPatchStart[jpatch];

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
		exit(1);
	}


	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType&, const VectorType&) const
	{
		err("BatchedGemm: matrixVector not implemented yet, sorry\n");
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
};
}
#endif // BATCHEDGEMM_H
