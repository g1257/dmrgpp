#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include <cassert>
#include <complex>
#include "Vector.h"
#include "../../../../dmrgppPluginSc/src/BatchedGemm.h"

typedef PsimagLite::Vector<IntegerType>::Type VectorIntegerType;
typedef PsimagLite::Vector<int>::Type VectorIntType;

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm2 {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename MatrixDenseOrSparseType::MatrixType MatrixType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename GenIjPatchType::BasisType BasisType;
	typedef BatchedGemm<ComplexOrRealType> BatchedGemmPluginScType;

	static const typename InitKronType::WhatBasisEnum DUMMY = InitKronType::OLD;

public:

	BatchedGemm2(const InitKronType& initKron)
	    :  initKron_(initKron),
	      Abatch_(0),
	      Bbatch_(0),
	      leftPatchStart_(0),
	      rightPatchStart_(0),
	      xyPatchStart_(0),
	      batchedGemm_(0)
	{
		if (!enabled()) return;
		convertOffsets(offsets_);
		SizeType npatches = initKron_.numberOfPatches(DUMMY);
		SizeType nC = initKron_.connections();
		ComplexOrRealType** aptr = new ComplexOrRealType*[npatches*npatches*nC];
		ComplexOrRealType** bptr = new ComplexOrRealType*[npatches*npatches*nC];
		VectorIntType ldAptr(npatches*npatches*nC);
		VectorIntType ldBptr(npatches*npatches*nC);

		convertToVector(pLeft_, initKron_.lrs(DUMMY).left());
		convertToVector(pRight_, initKron_.lrs(DUMMY).right());
		assert(pLeft_.size() == npatches);
		assert(pRight_.size() == npatches);

		for (SizeType outPatch = 0; outPatch < npatches; ++outPatch) {
			for (SizeType inPatch = 0; inPatch < npatches; ++inPatch) {
				for (SizeType ic=0;ic<nC;++ic) {
					const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
					const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

					const MatrixDenseOrSparseType& Amat = xiStruct(outPatch,inPatch);
					const MatrixDenseOrSparseType& Bmat = yiStruct(outPatch,inPatch);

					const MatrixType& AmatDense = Amat.dense();
					const MatrixType& BmatDense = Bmat.dense();

					ComplexOrRealType* a = const_cast<ComplexOrRealType*>(&(AmatDense(0,0)));
					ComplexOrRealType* b = const_cast<ComplexOrRealType*>(&(BmatDense(0,0)));
					aptr[outPatch + inPatch*npatches + ic*npatches*npatches] = a;
					bptr[outPatch + inPatch*npatches + ic*npatches*npatches] = b;

#ifndef NDEBUG
					SizeType nrowA = pLeft_[outPatch];
					SizeType ncolA = pLeft_[inPatch];
					SizeType nrowB = pRight_[outPatch];
					SizeType ncolB = pRight_[inPatch];

					assert( nrowA == AmatDense.rows() );
					assert( ncolA == AmatDense.cols() );
					assert( nrowB == BmatDense.rows() );
					assert( ncolB == BmatDense.cols() );
#endif

					ldAptr[outPatch + inPatch*npatches + ic*npatches*npatches] = AmatDense.rows();
					ldBptr[outPatch + inPatch*npatches + ic*npatches*npatches] = BmatDense.rows();
				}
			}
		}

		batchedGemm_ = new BatchedGemmPluginScType(nC,
		                                           npatches,
		                                           &(pLeft_[0]),
		        &(pRight_[0]),
		        aptr,
		        &(ldAptr[0]),
		        bptr,
		        &(ldBptr[0]));
		delete[] aptr;
		aptr = 0;
		delete[] bptr;
		bptr = 0;
	}

	BatchedGemm2()
	{
		delete batchedGemm_;
		batchedGemm_ = 0;

	}

	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		assert(enabled());
		ComplexOrRealType* vinptr = const_cast<ComplexOrRealType*>(&(vin[0]));
		ComplexOrRealType* voutptr = const_cast<ComplexOrRealType*>(&(vout[0]));
		batchedGemm_->apply_Htarget(vinptr, voutptr);
	}

private:

	void convertToVector(VectorIntType& v, const BasisType& b) const
	{
		SizeType total = b.partition() - 1;
		v.clear();
		v.resize(total, 0);
		for (SizeType i = 0; i < total; ++i)
			v[i] = b.partition(i + 1) - b.partition(i);
	}

	void convertOffsets(VectorIntegerType& v) const
	{
		SizeType npatches = initKron_.numberOfPatches(DUMMY);
		v.clear();
		v.resize(npatches, 0);
		for (SizeType i = 0; i < npatches; ++i)
			v[i] =initKron_.offsetForPatches(DUMMY, i);
	}

	BatchedGemm2(const BatchedGemm2&);

	BatchedGemm2& operator=(const BatchedGemm2&);

	const InitKronType& initKron_;
	VectorIntegerType offsets_;
	VectorIntType pLeft_;
	VectorIntType pRight_;
	mutable ComplexOrRealType* Abatch_;
	mutable ComplexOrRealType* Bbatch_;
	mutable IntegerType* leftPatchStart_;
	mutable IntegerType* rightPatchStart_;
	mutable IntegerType* xyPatchStart_;
	mutable int* ldAbatch_;
	mutable int* ldBbatch_;
	BatchedGemm<ComplexOrRealType>* batchedGemm_;
};
}
#endif // BATCHEDGEMM_H
