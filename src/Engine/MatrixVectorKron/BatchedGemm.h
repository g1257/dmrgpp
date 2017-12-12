#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include <complex>
#include "Vector.h"

typedef long int BatchedDgemmIntegerType;
typedef PsimagLite::Vector<BatchedDgemmIntegerType>::Type VectorBatchedDgemmIntegerType;
typedef double BatchedDgemmFpType;
typedef std::complex<BatchedDgemmFpType> BatchedDgemmComplexType;

/*** BEGIN Functions symbols that need to be provided by PLUGIN_SC */

// All OUTPUTS are for library internal use
extern "C" void setup_vbatch(BatchedDgemmIntegerType, // number of connections (INPUT)
                             BatchedDgemmIntegerType, // number of patches (INPUT)
                             BatchedDgemmIntegerType*, // left_patch_size[] (INPUT)
                             BatchedDgemmIntegerType*, // right_patch_size[] (INPUT)
                             BatchedDgemmIntegerType**, // left_patch_start[] (OUTPUT)
                             BatchedDgemmIntegerType**, // right_patch_start[] (OUTPUT)
                             BatchedDgemmIntegerType**, // xy_patch_start[] (OUTPUT)
                             BatchedDgemmFpType**, // Abatch[] (OUTPUT)
                             BatchedDgemmIntegerType**, // ld_Abatch Matrices A (OUTPUT)
                             BatchedDgemmFpType**, // Bbatch[] (OUTPUT)
                             BatchedDgemmIntegerType**, // ld_Bbatch (OUTPUT)
                             BatchedDgemmFpType**, // Matrices A (INPUT)
                             BatchedDgemmIntegerType*, // Rows of Matrices A (INPUT)
                             BatchedDgemmFpType**, // Matrices B (INPUT)
                             BatchedDgemmIntegerType*); // Rows of Matrices B (INPUT)

// All INPUTS except for vout, Abatch and Bbatch are for library internal use
extern "C" void apply_Htarget_vbatch(BatchedDgemmIntegerType, // number of connections
                                     BatchedDgemmIntegerType, // number of patches
                                     BatchedDgemmIntegerType*, // permutation left
                                     BatchedDgemmIntegerType*, // permutation right
                                     BatchedDgemmIntegerType*, // summed offsets
                                     BatchedDgemmFpType*, // Abatch for library internal use
                                     BatchedDgemmFpType*, // Bbatch for library internal use
                                     BatchedDgemmFpType*, // vin, already copied in
                                     BatchedDgemmFpType*); // vout, will be copied out (OUTPUT)

// SEE setup_vbatch above
// All OUTPUTS are for library internal use
extern "C" void setup_vbatchCmplx(BatchedDgemmIntegerType, // number of connections (INPUT)
                                  BatchedDgemmIntegerType, // number of patches (INPUT)
                                  BatchedDgemmIntegerType*, // left_patch_size[] (INPUT)
                                  BatchedDgemmIntegerType*, // right_patch_size[] (INPUT)
                                  BatchedDgemmIntegerType**, // left_patch_start[] (OUTPUT)
                                  BatchedDgemmIntegerType**, // right_patch_start[] (OUTPUT)
                                  BatchedDgemmIntegerType**, // xy_patch_start[] (OUTPUT)
                                  BatchedDgemmComplexType**, // Abatch[] (OUTPUT)
                                  BatchedDgemmIntegerType**, // ld_Abatch Matrices A (OUTPUT)
                                  BatchedDgemmComplexType**, // Bbatch[] (OUTPUT)
                                  BatchedDgemmIntegerType**, // ld_Bbatch (OUTPUT)
                                  BatchedDgemmComplexType**, // Matrices A (INPUT)
                                  BatchedDgemmIntegerType*, // Rows of Matrices A (INPUT)
                                  BatchedDgemmComplexType**, // Matrices B (INPUT)
                                  BatchedDgemmIntegerType*); // Rows of Matrices B (INPUT)

// SEE apply_Htarget_vbatch above
extern "C" void apply_Htarget_vbatchCmplx(BatchedDgemmIntegerType,
                                          BatchedDgemmIntegerType,
                                          BatchedDgemmIntegerType*,
                                          BatchedDgemmIntegerType*,
                                          BatchedDgemmIntegerType*,
                                          BatchedDgemmComplexType*,
                                          BatchedDgemmComplexType*,
                                          BatchedDgemmComplexType*,
                                          BatchedDgemmComplexType*);

// All OUTPUTS, all for library internal use
extern "C" void unsetup_vbatch(BatchedDgemmIntegerType**, // left_patch_start[] (OUTPUT)
                               BatchedDgemmIntegerType**, // right_patch_start[] (OUTPUT)
                               BatchedDgemmIntegerType**, // xy_patch_start[] (OUTPUT)
                               BatchedDgemmFpType**, // Abatch[] (OUTPUT)
                               BatchedDgemmIntegerType**, // ld_Abatch Matrices A (OUTPUT)
                               BatchedDgemmFpType**, // Bbatch[] (OUTPUT)
                               BatchedDgemmIntegerType**); // ld_Bbatch (OUTPUT)

// SEE unsetup_vbatch above
// All OUTPUTS, all for library internal use
extern "C" void unsetup_vbatchCmplx(BatchedDgemmIntegerType**, // left_patch_start[] (OUTPUT)
                                    BatchedDgemmIntegerType**, // right_patch_start[] (OUTPUT)
                                    BatchedDgemmIntegerType**, // xy_patch_start[] (OUTPUT)
                                    BatchedDgemmComplexType**, // Abatch[] (OUTPUT)
                                    BatchedDgemmIntegerType**, // ld_Abatch Matrices A (OUTPUT)
                                    BatchedDgemmComplexType**, // Bbatch[] (OUTPUT)
                                    BatchedDgemmIntegerType**); // ld_Bbatch (OUTPUT)

/*** END Functions symbols that need to be provided by PLUGIN_SC */

template<typename T>
void applyHtargetVbatch(BatchedDgemmIntegerType,
                        BatchedDgemmIntegerType,
                        const VectorBatchedDgemmIntegerType&,
                        const VectorBatchedDgemmIntegerType&,
                        const VectorBatchedDgemmIntegerType&,
                        T*,
                        T*,
                        T*,
                        T*);

template<>
inline void applyHtargetVbatch<BatchedDgemmComplexType>(BatchedDgemmIntegerType a,
                                                        BatchedDgemmIntegerType b,
                                                        const VectorBatchedDgemmIntegerType& c,
                                                        const VectorBatchedDgemmIntegerType& d,
                                                        const VectorBatchedDgemmIntegerType& e,
                                                        BatchedDgemmComplexType* f,
                                                        BatchedDgemmComplexType* g,
                                                        BatchedDgemmComplexType* h,
                                                        BatchedDgemmComplexType* i)
{
	BatchedDgemmIntegerType* cptr = const_cast<BatchedDgemmIntegerType*>(&(c[0]));
	BatchedDgemmIntegerType* dptr = const_cast<BatchedDgemmIntegerType*>(&(d[0]));
	BatchedDgemmIntegerType* eptr = const_cast<BatchedDgemmIntegerType*>(&(e[0]));
	apply_Htarget_vbatchCmplx(a, b, cptr, dptr, eptr, f, g, h, i);
}

template<>
inline void applyHtargetVbatch<BatchedDgemmFpType>(BatchedDgemmIntegerType a,
                                                   BatchedDgemmIntegerType b,
                                                   const VectorBatchedDgemmIntegerType& c,
                                                   const VectorBatchedDgemmIntegerType& d,
                                                   const VectorBatchedDgemmIntegerType& e,
                                                   BatchedDgemmFpType* f,
                                                   BatchedDgemmFpType* g,
                                                   BatchedDgemmFpType* h,
                                                   BatchedDgemmFpType* i)
{
	BatchedDgemmIntegerType* cptr = const_cast<BatchedDgemmIntegerType*>(&(c[0]));
	BatchedDgemmIntegerType* dptr = const_cast<BatchedDgemmIntegerType*>(&(d[0]));
	BatchedDgemmIntegerType* eptr = const_cast<BatchedDgemmIntegerType*>(&(e[0]));
	apply_Htarget_vbatch(a, b, cptr, dptr, eptr, f, g, h, i);
}

/******/

template<typename T>
void setupVbatch(BatchedDgemmIntegerType,
                 BatchedDgemmIntegerType,
                 const VectorBatchedDgemmIntegerType&,
                 const VectorBatchedDgemmIntegerType&,
                 BatchedDgemmIntegerType**,
                 BatchedDgemmIntegerType**,
                 BatchedDgemmIntegerType**,
                 T**,
                 BatchedDgemmIntegerType**,
                 T**,
                 BatchedDgemmIntegerType**,
                 T**,
                 const VectorBatchedDgemmIntegerType&,
                 T**,
                 const VectorBatchedDgemmIntegerType&);

template<>
inline void setupVbatch<BatchedDgemmComplexType>(BatchedDgemmIntegerType a,
                                                 BatchedDgemmIntegerType b,
                                                 const VectorBatchedDgemmIntegerType& c,
                                                 const VectorBatchedDgemmIntegerType& d,
                                                 BatchedDgemmIntegerType** e,
                                                 BatchedDgemmIntegerType** f,
                                                 BatchedDgemmIntegerType** g,
                                                 BatchedDgemmComplexType** h,
                                                 BatchedDgemmIntegerType** i,
                                                 BatchedDgemmComplexType** j,
                                                 BatchedDgemmIntegerType** k,
                                                 BatchedDgemmComplexType** l,
                                                 const VectorBatchedDgemmIntegerType& m,
                                                 BatchedDgemmComplexType** n,
                                                 const VectorBatchedDgemmIntegerType& o)
{
	BatchedDgemmIntegerType* cptr = const_cast<BatchedDgemmIntegerType*>(&(c[0]));
	BatchedDgemmIntegerType* dptr = const_cast<BatchedDgemmIntegerType*>(&(d[0]));
	BatchedDgemmIntegerType* mptr = const_cast<BatchedDgemmIntegerType*>(&(m[0]));
	BatchedDgemmIntegerType* optr = const_cast<BatchedDgemmIntegerType*>(&(o[0]));
	setup_vbatchCmplx(a, b, cptr, dptr, e, f, g, h, i, j, k, l, mptr, n, optr);
}

template<>
inline void setupVbatch<BatchedDgemmFpType>(BatchedDgemmIntegerType a,
                                            BatchedDgemmIntegerType b,
                                            const VectorBatchedDgemmIntegerType& c,
                                            const VectorBatchedDgemmIntegerType& d,
                                            BatchedDgemmIntegerType** e,
                                            BatchedDgemmIntegerType** f,
                                            BatchedDgemmIntegerType** g,
                                            BatchedDgemmFpType** h,
                                            BatchedDgemmIntegerType** i,
                                            BatchedDgemmFpType** j,
                                            BatchedDgemmIntegerType** k,
                                            BatchedDgemmFpType** l,
                                            const VectorBatchedDgemmIntegerType& m,
                                            BatchedDgemmFpType** n,
                                            const VectorBatchedDgemmIntegerType& o)
{
	BatchedDgemmIntegerType* cptr = const_cast<BatchedDgemmIntegerType*>(&(c[0]));
	BatchedDgemmIntegerType* dptr = const_cast<BatchedDgemmIntegerType*>(&(d[0]));
	BatchedDgemmIntegerType* mptr = const_cast<BatchedDgemmIntegerType*>(&(m[0]));
	BatchedDgemmIntegerType* optr = const_cast<BatchedDgemmIntegerType*>(&(o[0]));
	setup_vbatch(a, b, cptr, dptr, e, f, g, h, i, j, k, l, mptr, n, optr);
}

template<typename T>
void unsetupVbatch(BatchedDgemmIntegerType**,
                   BatchedDgemmIntegerType**,
                   BatchedDgemmIntegerType**,
                   T**,
                   BatchedDgemmIntegerType**,
                   T**,
                   BatchedDgemmIntegerType**);

template<>
inline void unsetupVbatch<BatchedDgemmFpType>(BatchedDgemmIntegerType** a,
                                              BatchedDgemmIntegerType** b,
                                              BatchedDgemmIntegerType** c,
                                              BatchedDgemmFpType** d,
                                              BatchedDgemmIntegerType** e,
                                              BatchedDgemmFpType** f,
                                              BatchedDgemmIntegerType** g)
{
	unsetup_vbatch(a, b, c, d, e, f, g);
}

template<>
inline void unsetupVbatch<BatchedDgemmComplexType>(BatchedDgemmIntegerType** a,
                                                   BatchedDgemmIntegerType** b,
                                                   BatchedDgemmIntegerType** c,
                                                   BatchedDgemmComplexType** d,
                                                   BatchedDgemmIntegerType** e,
                                                   BatchedDgemmComplexType** f,
                                                   BatchedDgemmIntegerType** g)
{
	unsetup_vbatchCmplx(a, b, c, d, e, f, g);
}

/******/

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef typename InitKronType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename MatrixDenseOrSparseType::MatrixType MatrixType;
	typedef typename InitKronType::GenIjPatchType GenIjPatchType;
	typedef typename GenIjPatchType::BasisType BasisType;

	// Say 1 for 1-based, 0 for 0-based
	static const SizeType baseForIntegerVectors_ = 1;

public:

	BatchedGemm(const InitKronType& initKron)
	    :  initKron_(initKron),
	      Abatch_(0),
	      Bbatch_(0),
	      leftPatchStart_(0),
	      rightPatchStart_(0),
	      xyPatchStart_(0)
	{
		convertOffsets(offsets_);
		SizeType total = initKron_.numberOfPatches(InitKronType::OLD);
		SizeType nC = initKron_.connections();
		ComplexOrRealType** aptr = new ComplexOrRealType*[total*total*nC];
		ComplexOrRealType** bptr = new ComplexOrRealType*[total*total*nC];
		VectorBatchedDgemmIntegerType ldAptr(total*total*nC);
		VectorBatchedDgemmIntegerType ldBptr(total*total*nC);

		for (SizeType outPatch = 0; outPatch < total; ++outPatch) {
			for (SizeType inPatch = 0; inPatch < total; ++inPatch) {
				for (SizeType ic=0;ic<nC;++ic) {
					const ArrayOfMatStructType& xiStruct = initKron_.xc(ic);
					const ArrayOfMatStructType& yiStruct = initKron_.yc(ic);

					const MatrixDenseOrSparseType& Amat =  xiStruct(outPatch,inPatch);
					const MatrixDenseOrSparseType& Bmat =  yiStruct(outPatch,inPatch);

					const MatrixType& AmatDense = Amat.dense();
					const MatrixType& BmatDense = Bmat.dense();

					ComplexOrRealType* a = const_cast<ComplexOrRealType*>(&(AmatDense(0,0)));
					ComplexOrRealType* b = const_cast<ComplexOrRealType*>(&(BmatDense(0,0)));
					aptr[inPatch + outPatch*total + ic*total*total] = a;
					bptr[inPatch + outPatch*total + ic*total*total] = b;

					ldAptr[inPatch + outPatch*total + ic*total*total] = AmatDense.rows();
					ldBptr[inPatch + outPatch*total + ic*total*total] = BmatDense.rows();
				}
			}
		}

		// call setup_vbatch and fill aBatch_ and bBatch_

		convertToVector(pLeft_, initKron_.lrs(InitKronType::NEW).left());
		convertToVector(pRight_, initKron_.lrs(InitKronType::NEW).right());

		setupVbatch(initKron.connections(), // number of connections (INPUT)
		            initKron_.patch(InitKronType::NEW, GenIjPatchType::LEFT).size(),
		            pLeft_,
		            pRight_,
		            &leftPatchStart_,
		            &rightPatchStart_,
		            &xyPatchStart_,
		            &Abatch_,
		            &ldAbatch_,
		            &Bbatch_,
		            &ldBbatch_,
		            aptr,
		            ldAptr,
		            bptr,
		            ldBptr);

		delete[] aptr;
		aptr = 0;
		delete[] bptr;
		bptr = 0;
	}

	~BatchedGemm()
	{
		unsetupVbatch(&leftPatchStart_,
		              &rightPatchStart_,
		              &xyPatchStart_,
		              &Abatch_,
		              &ldAbatch_,
		              &Bbatch_,
		              &ldBbatch_);
	}

	bool enabled() const { return true; }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		ComplexOrRealType* vinptr = const_cast<ComplexOrRealType*>(&(vin[0]));
		ComplexOrRealType* voutptr = const_cast<ComplexOrRealType*>(&(vout[0]));
		applyHtargetVbatch(initKron_.connections(),
		                   initKron_.patch(InitKronType::NEW, GenIjPatchType::LEFT).size(),
		                   pLeft_,
		                   pRight_,
		                   offsets_,
		                   Abatch_,
		                   Bbatch_,
		                   vinptr,
		                   voutptr);
	}

private:

	void convertToVector(VectorBatchedDgemmIntegerType& v, const BasisType& b) const
	{
		SizeType total = b.partition();
		v.clear();
		v.resize(total + baseForIntegerVectors_, 0);
		for (SizeType i = 0; i < total; ++i)
			v[i + baseForIntegerVectors_] = b.partition(i);
	}

	void convertOffsets(VectorBatchedDgemmIntegerType& v) const
	{
		SizeType total = initKron_.offsetForPatches(InitKronType::NEW);
		v.clear();
		v.resize(total + baseForIntegerVectors_, 0);
		for (SizeType i = 0; i < total; ++i)
			v[i + baseForIntegerVectors_] =initKron_.offsetForPatches(InitKronType::NEW, i);
	}

	BatchedGemm(const BatchedGemm&);

	BatchedGemm& operator=(const BatchedGemm&);

	const InitKronType& initKron_;
	VectorBatchedDgemmIntegerType offsets_;
	VectorBatchedDgemmIntegerType pLeft_;
	VectorBatchedDgemmIntegerType pRight_;
	ComplexOrRealType* Abatch_;
	ComplexOrRealType* Bbatch_;
	BatchedDgemmIntegerType* leftPatchStart_;
	BatchedDgemmIntegerType* rightPatchStart_;
	BatchedDgemmIntegerType* xyPatchStart_;
	BatchedDgemmIntegerType* ldAbatch_;
	BatchedDgemmIntegerType* ldBbatch_;
};
}
#endif // BATCHEDGEMM_H
