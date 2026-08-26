#ifndef DMRGPPPLUGINSC_H
#define DMRGPPPLUGINSC_H

#include "dmrg_types.h"
#include <PsimagLite/Matrix.h>
#include <PsimagLite/Vector.h>

#include <Kokkos_Core.hpp>
#include <PsimagLite/KokkosType.h>

template <typename T> class BatchedGemm {

	typedef int            IntegerType;
	typedef std::vector<T> VectorType;
	using KokkosScalar = PsimagLite::KokkosType<T>::type;
	typedef PsimagLite::Vector<SizeType>::Type    VectorSizeType;
	typedef PsimagLite::Vector<IntegerType>::Type VectorIntegerType;

public:

	BatchedGemm(SizeType         noperator,
	            SizeType         npatches,
	            VectorSizeType&  left_patch_size,
	            VectorSizeType&  right_patch_size,
	            std::vector<T*>& Amatrix,
	            VectorSizeType&  ld_Amatrix,
	            std::vector<T*>& Bmatrix,
	            VectorSizeType&  ld_Bmatrix);

	void apply_Htarget(Kokkos::View<const KokkosScalar*, Kokkos::SharedSpace> vin,
	                   Kokkos::View<KokkosScalar*, Kokkos::SharedSpace>       vout);

	~BatchedGemm();

private:

	BatchedGemm(const BatchedGemm<T>&);
	BatchedGemm& operator=(const BatchedGemm<T>&);

	SizeType noperator_;
	SizeType npatches_;
	SizeType ld_Abatch_;
	SizeType ld_Bbatch_;

	VectorType Abatch_;
	VectorType Bbatch_;
	// Set up in setup_batch.cpp
	VectorSizeType    left_patch_start_;
	VectorSizeType    right_patch_start_;
	VectorSizeType    xy_patch_start_;
	VectorSizeType    nC_;
	VectorIntegerType ld_gAbatch_;
	VectorIntegerType ld_gBbatch_;
	T**               gAbatch_;
	T**               gBbatch_;
};

#endif
