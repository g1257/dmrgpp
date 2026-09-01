#ifndef BATCHED_GEMM_PLUGIN_SC_H
#define BATCHED_GEMM_PLUGIN_SC_H
// Don't include this file directly; use BatchedGemmInclude.hh

#include "BatchedGemm.h"
#include <PsimagLite/KokkosType.h>
#include <PsimagLite/Matrix.h>
#include <PsimagLite/ProgressIndicator.h>
#include <PsimagLite/Vector.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

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
		setup_();
	}

	void setup_()
	{
		Kokkos::Profiling::ScopedRegion region("BatchedGemmPluginSc::setup");

		SizeType                        npatches = initKron_.numberOfPatches(DUMMY);
		SizeType                        nC       = initKron_.connections();
		const SizeType                  total    = npatches * npatches * nC;
		std::vector<ComplexOrRealType*> aptr(total, 0);
		std::vector<ComplexOrRealType*> bptr(total, 0);
		VectorSizeType                  ldAptr(npatches * npatches * nC);
		VectorSizeType                  ldBptr(npatches * npatches * nC);

		pLeft_.resize(npatches, 0);
		pRight_.resize(npatches, 0);

		{
			Kokkos::Profiling::ScopedRegion region(
			    "BatchedGemmPluginSc::setup::patch_pointers");
			SizeType zeroes = 0;
			for (SizeType ic = 0; ic < nC; ++ic) {
				for (SizeType inPatch = 0; inPatch < npatches; ++inPatch) {
					for (SizeType outPatch = 0; outPatch < npatches;
					     ++outPatch) {

						const ArrayOfMatStructType& xiStruct
						    = initKron_.xc(ic);
						const ArrayOfMatStructType& yiStruct
						    = initKron_.yc(ic);

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
		}

		PsimagLite::OstringStream msg(std::cout.precision());
		msg() << "PLUGIN_SC: is in use, npatches=" << npatches;
		msg() << " connections=" << nC << " zeroConnections=" << zeroes;
		progress_.printline(msg, std::cout);

		batchedGemm_ = new BatchedGemmType(
		    nC, npatches, pLeft_, pRight_, aptr, ldAptr, bptr, ldBptr);
	}

	~BatchedGemmPluginSc()
	{
		delete batchedGemm_;
		batchedGemm_ = nullptr;
		for (SizeType i = 0; i < garbage_.size(); ++i) {
			delete garbage_[i];
			garbage_[i] = nullptr;
		}
	}

	bool enabled() const { return initKron_.params().options.isSet("BatchedGemm"); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		Kokkos::Profiling::ScopedRegion region("BatchedGemmPluginSc::matrixVector");

		assert(enabled());
		using KokkosScalar = PsimagLite::KokkosType<ComplexOrRealType>::type;
		Kokkos::DefaultExecutionSpace exec;
		Kokkos::SharedSpace           mem;

		Kokkos::View<const KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
		     vin_view(reinterpret_cast<const KokkosScalar*>(vin.data()), vin.size());
		auto vin_view_shared
		    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), vin_view);

		Kokkos::View<KokkosScalar*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> vout_view(
		    reinterpret_cast<KokkosScalar*>(vout.data()), vout.size());
		auto vout_view_shared
		    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), vout_view);
		exec.fence();

		batchedGemm_->apply_Htarget(vin_view_shared, vout_view_shared);
		Kokkos::deep_copy(vout_view, vout_view_shared);
	}

private:

	void getMatrixPointers(ComplexOrRealType**            a,
	                       ComplexOrRealType**            b,
	                       const MatrixDenseOrSparseType& Amat,
	                       const MatrixDenseOrSparseType& Bmat) const
	{
		*a = nullptr;
		*b = nullptr;
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
