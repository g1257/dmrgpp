#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include <cassert>
#include <complex>
#include "Matrix.h"
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
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Vector<MatrixType*>::Type VectorMatrixType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename GenIjPatchType::BasisType BasisType;
	typedef BatchedGemm<ComplexOrRealType> BatchedGemmPluginScType;

	static const typename InitKronType::WhatBasisEnum DUMMY = InitKronType::OLD;

public:

	BatchedGemm2(const InitKronType& initKron)
	    :  progress_("BatchedGemm"), initKron_(initKron), batchedGemm_(0)
	{
		if (!enabled()) return;
		SizeType npatches = initKron_.numberOfPatches(DUMMY);
		SizeType nC = initKron_.connections();
		const SizeType total = npatches*npatches*nC;
		ComplexOrRealType** aptr = new ComplexOrRealType*[total];
		ComplexOrRealType** bptr = new ComplexOrRealType*[total];
		VectorIntType ldAptr(npatches*npatches*nC);
		VectorIntType ldBptr(npatches*npatches*nC);

		const ComplexOrRealType* zero = 0;
		memcpy(aptr, &zero, total*sizeof(ComplexOrRealType*));
		memcpy(bptr, &zero, total*sizeof(ComplexOrRealType*));

		pLeft_.resize(npatches, 0);
		pRight_.resize(npatches, 0);

		SizeType zeroes = 0;
		for (SizeType ic = 0; ic < nC; ++ic) {
			for (SizeType inPatch = 0; inPatch < npatches; ++inPatch) {
				for (SizeType outPatch = 0; outPatch < npatches; ++outPatch) {

					const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
					const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

					const MatrixDenseOrSparseType* Amat = xiStruct(outPatch,inPatch);
					const MatrixDenseOrSparseType* Bmat = yiStruct(outPatch,inPatch);

					if (!Amat || !Bmat) continue;

					ComplexOrRealType* a = 0;
					ComplexOrRealType* b = 0;
					getMatrixPointers(&a, &b, *Amat, *Bmat);

					if (a == 0) {
						assert(b == 0);
						++zeroes;
					}

					aptr[outPatch + inPatch*npatches + ic*npatches*npatches] = a;
					bptr[outPatch + inPatch*npatches + ic*npatches*npatches] = b;

					initKron_.checks(*Amat, *Bmat, outPatch, inPatch);
					pLeft_[inPatch] = Amat->cols();
					pRight_[inPatch] = Bmat->cols();

					ldAptr[outPatch + inPatch*npatches + ic*npatches*npatches] = Amat->rows();
					ldBptr[outPatch + inPatch*npatches + ic*npatches*npatches] = Bmat->rows();
				}
			}
		}

		{
			PsimagLite::OstringStream msg;
			msg<<"PLUGIN_SC: is in use, npatches="<<npatches;
			msg<<" connections="<<nC<<" zeroConnections="<<zeroes;
			progress_.printline(msg,std::cout);
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

	~BatchedGemm2()
	{
		delete batchedGemm_;
		batchedGemm_ = 0;
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
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

	void getMatrixPointers(ComplexOrRealType** a,
	                       ComplexOrRealType** b,
	                       const MatrixDenseOrSparseType& Amat,
	                       const MatrixDenseOrSparseType& Bmat) const
	{
		*a = *b = 0;
		if (Amat.isZero() || Bmat.isZero()) return;

		*a = getMatrixPointer(Amat);
		*b = getMatrixPointer(Bmat);
	}

	ComplexOrRealType* getMatrixPointer(const MatrixDenseOrSparseType& mat) const
	{
		if (!mat.isDense()) {
			MatrixType* matDense = new MatrixType();
			crsMatrixToFullMatrix(*matDense, mat.sparse());
			garbage_.push_back(matDense);
			return const_cast<ComplexOrRealType*>(&(matDense->operator()(0,0)));
		}

		return const_cast<ComplexOrRealType*>(&(mat.dense()(0,0)));
	}

	BatchedGemm2(const BatchedGemm2&);

	BatchedGemm2& operator=(const BatchedGemm2&);

	PsimagLite::ProgressIndicator progress_;
	const InitKronType& initKron_;
	VectorIntType pLeft_;
	VectorIntType pRight_;
	BatchedGemm<ComplexOrRealType>* batchedGemm_;
	mutable VectorMatrixType garbage_;
};
}
#endif // BATCHEDGEMM_H
