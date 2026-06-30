#ifndef BATCHED_GEMM_PLUGIN_SC_H
#define BATCHED_GEMM_PLUGIN_SC_H
// Don't include this file directly; use BatchedGemmInclude.hh

#include "BatchedGemm.h"
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "Vector.h"
#include <cassert>
#include <complex>

using VectorIntegerType = PsimagLite::Vector<IntegerType>::Type;
using VectorIntType     = PsimagLite::Vector<int>::Type;

namespace Dmrg {

template <typename InitKronType> class BatchedGemmPluginSc {

	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using VectorType              = typename MatrixDenseOrSparseType::VectorType;
	using VectorSizeType          = PsimagLite::Vector<SizeType>::Type;
	using SparseMatrixType        = typename InitKronType::SparseMatrixType;
	using ComplexOrRealType       = typename SparseMatrixType::value_type;
	using MatrixType              = PsimagLite::Matrix<ComplexOrRealType>;
	using VectorMatrixType        = typename PsimagLite::Vector<MatrixType*>::Type;
	using GenIjPatchType          = typename InitKronType::GenIjPatchType;
	using BasisType               = typename GenIjPatchType::BasisType;
	using BatchedGemmType         = BatchedGemm<ComplexOrRealType>;

	static const typename InitKronType::WhatBasisEnum DUMMY = InitKronType::OLD;

public:

	BatchedGemmPluginSc(const InitKronType& initKron)
	    : progress_("BatchedGemm")
	    , initKron_(initKron)
	    , batchedGemm_(0)
	{
		if (!enabled())
			return;
		SizeType                        npatches = initKron_.numberOfPatches(DUMMY);
		SizeType                        nC       = initKron_.connections();
		const SizeType                  total    = npatches * npatches * nC;
		std::vector<ComplexOrRealType*> aptr(total, 0);
		std::vector<ComplexOrRealType*> bptr(total, 0);
		VectorSizeType                  ldAptr(npatches * npatches * nC);
		VectorSizeType                  ldBptr(npatches * npatches * nC);

		pLeft_.resize(npatches, 0);
		pRight_.resize(npatches, 0);

		SizeType zeroes = 0;
		for (SizeType ic = 0; ic < nC; ++ic) {
			for (SizeType inPatch = 0; inPatch < npatches; ++inPatch) {
				for (SizeType outPatch = 0; outPatch < npatches; ++outPatch) {

					const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
					const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

					const MatrixDenseOrSparseType* Amat
					    = xiStruct(outPatch, inPatch);
					const MatrixDenseOrSparseType* Bmat
					    = yiStruct(outPatch, inPatch);

					if (!Amat || !Bmat)
						continue;

					ComplexOrRealType* a = 0;
					ComplexOrRealType* b = 0;
					getMatrixPointers(&a, &b, *Amat, *Bmat);

					if (a == 0) {
						assert(b == 0);
						++zeroes;
					}

					aptr[outPatch + inPatch * npatches
					     + ic * npatches * npatches]
					    = a;
					bptr[outPatch + inPatch * npatches
					     + ic * npatches * npatches]
					    = b;

					initKron_.checks(*Amat, *Bmat, outPatch, inPatch);
					pLeft_[inPatch]  = Amat->cols();
					pRight_[inPatch] = Bmat->cols();

					ldAptr[outPatch + inPatch * npatches
					       + ic * npatches * npatches]
					    = Amat->rows();
					ldBptr[outPatch + inPatch * npatches
					       + ic * npatches * npatches]
					    = Bmat->rows();
				}
			}
		}

		{
			PsimagLite::OstringStream msg(std::cout.precision());
			msg() << "PLUGIN_SC: is in use, npatches=" << npatches;
			msg() << " connections=" << nC << " zeroConnections=" << zeroes;
			progress_.printline(msg, std::cout);
		}

		batchedGemm_ = new BatchedGemmType(
		    nC, npatches, pLeft_, pRight_, aptr, ldAptr, bptr, ldBptr);
	}

	~BatchedGemmPluginSc()
	{
		delete batchedGemm_;
		batchedGemm_ = 0;
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			delete garbage_[i];
			garbage_[i] = 0;
		}
	}

	bool enabled() const { return initKron_.params().options.isSet("BatchedGemm"); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		assert(enabled());
		ComplexOrRealType* vinptr  = const_cast<ComplexOrRealType*>(&(vin[0]));
		ComplexOrRealType* voutptr = const_cast<ComplexOrRealType*>(&(vout[0]));
		batchedGemm_->apply_Htarget(vinptr, voutptr);
	}

private:

	void getMatrixPointers(ComplexOrRealType**            a,
	                       ComplexOrRealType**            b,
	                       const MatrixDenseOrSparseType& Amat,
	                       const MatrixDenseOrSparseType& Bmat) const
	{
		*a = *b = 0;
		if (Amat.isZero() || Bmat.isZero())
			return;

		*a = getMatrixPointer(Amat);
		*b = getMatrixPointer(Bmat);
	}

	ComplexOrRealType* getMatrixPointer(const MatrixDenseOrSparseType& mat) const
	{
		if (!mat.isDense()) {
			MatrixType* matDense = new MatrixType();
			crsMatrixToFullMatrix(*matDense, mat.sparse());
			garbage_.push_back(matDense);
			return const_cast<ComplexOrRealType*>(&(matDense->operator()(0, 0)));
		}

		return const_cast<ComplexOrRealType*>(&(mat.dense()(0, 0)));
	}

	BatchedGemmPluginSc(const BatchedGemmPluginSc&) = delete;

	BatchedGemmPluginSc& operator=(const BatchedGemmPluginSc&) = delete;

	PsimagLite::ProgressIndicator   progress_;
	const InitKronType&             initKron_;
	VectorSizeType                  pLeft_;
	VectorSizeType                  pRight_;
	BatchedGemm<ComplexOrRealType>* batchedGemm_;
	mutable VectorMatrixType        garbage_;
};
}
#endif // BATCHED_GEMM_PLUGIN_SC_H
