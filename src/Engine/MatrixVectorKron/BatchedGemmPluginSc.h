#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include <complex>
#include "Vector.h"
#include "../../../../dmrgppPluginSc/include/DmrgPluginSc.h"

//typedef long int IntegerType;
typedef PsimagLite::Vector<IntegerType>::Type VectorIntegerType;
//typedef double FpType;
typedef std::complex<FpType> BatchedDgemmComplexType;

/*** BEGIN Functions symbols that need to be provided by PLUGIN_SC */
#if 0
// All OUTPUTS are for library internal use
extern "C" void setup_vbatch(IntegerType, // number of connections (INPUT)
                             IntegerType, // number of patches (INPUT)
                             IntegerType*, // left_patch_size[] (INPUT)
                             IntegerType*, // right_patch_size[] (INPUT)
                             IntegerType**, // left_patch_start[] (OUTPUT)
                             IntegerType**, // right_patch_start[] (OUTPUT)
                             IntegerType**, // xy_patch_start[] (OUTPUT)
                             FpType**, // Abatch[] (OUTPUT)
                             IntegerType**, // ld_Abatch Matrices A (OUTPUT)
                             FpType**, // Bbatch[] (OUTPUT)
                             IntegerType**, // ld_Bbatch (OUTPUT)
                             FpType**, // Matrices A (INPUT)
                             IntegerType*, // Rows of Matrices A (INPUT)
                             FpType**, // Matrices B (INPUT)
                             IntegerType*); // Rows of Matrices B (INPUT)

// OUTPUT vout is used
// All other OUTPUTS are for library internal use
extern "C" void apply_Htarget_vbatch(IntegerType, // number of connections (INPUT)
                                     IntegerType, // number of patches (INPUT)
                                     IntegerType*, // left_patch_size[] (INPUT)
                                     IntegerType*, // right_patch_size[] (INPUT)
                                     IntegerType**, // left_patch_start[] (OUTPUT)
                                     IntegerType**, // right_patch_start[] (OUTPUT)
                                     IntegerType**, // xy_patch_start[] (OUTPUT)
                                     FpType**, // Abatch[] (OUTPUT)
                                     IntegerType**, // ld_Abatch Matrices A (OUTPUT)
                                     FpType**, // Bbatch[] (OUTPUT)
                                     IntegerType**, // ld_Bbatch (OUTPUT)
                                     FpType*, // vin, already copied in
                                     FpType*); // vout, will be copied out (OUTPUT)

// SEE setup_vbatch above
// All OUTPUTS are for library internal use
extern "C" void setup_vbatchCmplx(IntegerType, // number of connections (INPUT)
                                  IntegerType, // number of patches (INPUT)
                                  IntegerType*, // left_patch_size[] (INPUT)
                                  IntegerType*, // right_patch_size[] (INPUT)
                                  IntegerType**, // left_patch_start[] (OUTPUT)
                                  IntegerType**, // right_patch_start[] (OUTPUT)
                                  IntegerType**, // xy_patch_start[] (OUTPUT)
                                  BatchedDgemmComplexType**, // Abatch[] (OUTPUT)
                                  IntegerType**, // ld_Abatch Matrices A (OUTPUT)
                                  BatchedDgemmComplexType**, // Bbatch[] (OUTPUT)
                                  IntegerType**, // ld_Bbatch (OUTPUT)
                                  BatchedDgemmComplexType**, // Matrices A (INPUT)
                                  IntegerType*, // Rows of Matrices A (INPUT)
                                  BatchedDgemmComplexType**, // Matrices B (INPUT)
                                  IntegerType*); // Rows of Matrices B (INPUT)

// SEE apply_Htarget_vbatch above
// OUTPUT vout is used
// All other OUTPUTS are for library internal use
extern "C" void apply_Htarget_vbatchCmplx(IntegerType, // number of connections (INPUT)
                                          IntegerType, // number of patches (INPUT)
                                          IntegerType*, // left_patch_size[] (INPUT)
                                          IntegerType*, // right_patch_size[] (INPUT)
                                          IntegerType**, // left_patch_start[] (OUTPUT)
                                          IntegerType**, // right_patch_start[] (OUTPUT)
                                          IntegerType**, // xy_patch_start[] (OUTPUT)
                                          BatchedDgemmComplexType**, // Abatch[] (OUTPUT)
                                          IntegerType**, // ld_Abatch Matrices A (OUTPUT)
                                          BatchedDgemmComplexType**, // Bbatch[] (OUTPUT)
                                          IntegerType**, // ld_Bbatch (OUTPUT)
                                          BatchedDgemmComplexType*,
                                          BatchedDgemmComplexType*);

// All OUTPUTS, all for library internal use
extern "C" void unsetup_vbatch(IntegerType**, // left_patch_start[] (OUTPUT)
                               IntegerType**, // right_patch_start[] (OUTPUT)
                               IntegerType**, // xy_patch_start[] (OUTPUT)
                               FpType**, // Abatch[] (OUTPUT)
                               IntegerType**, // ld_Abatch Matrices A (OUTPUT)
                               FpType**, // Bbatch[] (OUTPUT)
                               IntegerType**); // ld_Bbatch (OUTPUT)

// SEE unsetup_vbatch above
// All OUTPUTS, all for library internal use
extern "C" void unsetup_vbatchCmplx(IntegerType**, // left_patch_start[] (OUTPUT)
                                    IntegerType**, // right_patch_start[] (OUTPUT)
                                    IntegerType**, // xy_patch_start[] (OUTPUT)
                                    BatchedDgemmComplexType**, // Abatch[] (OUTPUT)
                                    IntegerType**, // ld_Abatch Matrices A (OUTPUT)
                                    BatchedDgemmComplexType**, // Bbatch[] (OUTPUT)
                                    IntegerType**); // ld_Bbatch (OUTPUT)

/*** END Functions symbols that need to be provided by PLUGIN_SC */
#endif

template<typename T>
void applyHtargetVbatch(IntegerType,
                        IntegerType,
                        const VectorIntegerType&,
                        const VectorIntegerType&,
                        IntegerType**,
                        IntegerType**,
                        IntegerType**,
                        T**,
                        IntegerType**,
                        T**,
                        IntegerType**,
                        T*,
                        T*);

template<>
inline void applyHtargetVbatch<BatchedDgemmComplexType>(IntegerType a,
                                                        IntegerType b,
                                                        const VectorIntegerType& c,
                                                        const VectorIntegerType& d,
                                                        IntegerType** e,
                                                        IntegerType** f,
                                                        IntegerType** g,
                                                        BatchedDgemmComplexType** h,
                                                        IntegerType** i,
                                                        BatchedDgemmComplexType** j,
                                                        IntegerType** k,
                                                        BatchedDgemmComplexType* l,
                                                        BatchedDgemmComplexType* m)
{
	//IntegerType* cptr = const_cast<IntegerType*>(&(c[0]));
	//IntegerType* dptr = const_cast<IntegerType*>(&(d[0]));
	err("PluginSc does not support complex\n");
	//apply_Htarget_vbatchCmplx(a, b, cptr, dptr, e, f, g, h, i, j, k, l, m);
}

template<>
inline void applyHtargetVbatch<FpType>(IntegerType a,
                                       IntegerType b,
                                       const VectorIntegerType& c,
                                       const VectorIntegerType& d,
                                       IntegerType** e,
                                       IntegerType** f,
                                       IntegerType** g,
                                       FpType** h,
                                       IntegerType** i,
                                       FpType** j,
                                       IntegerType** k,
                                       FpType* l,
                                       FpType* m)
{
	IntegerType* cptr = const_cast<IntegerType*>(&(c[0]));
	IntegerType* dptr = const_cast<IntegerType*>(&(d[0]));
	apply_Htarget_vbatch(a, b, cptr, dptr, e, f, g, h, i, j, k, l, m);
}

/******/

template<typename T>
void setupVbatch(IntegerType,
                 IntegerType,
                 const VectorIntegerType&,
                 const VectorIntegerType&,
                 IntegerType**,
                 IntegerType**,
                 IntegerType**,
                 T**,
                 IntegerType**,
                 T**,
                 IntegerType**,
                 T**,
                 const VectorIntegerType&,
                 T**,
                 const VectorIntegerType&);

template<>
inline void setupVbatch<BatchedDgemmComplexType>(IntegerType a,
                                                 IntegerType b,
                                                 const VectorIntegerType& c,
                                                 const VectorIntegerType& d,
                                                 IntegerType** e,
                                                 IntegerType** f,
                                                 IntegerType** g,
                                                 BatchedDgemmComplexType** h,
                                                 IntegerType** i,
                                                 BatchedDgemmComplexType** j,
                                                 IntegerType** k,
                                                 BatchedDgemmComplexType** l,
                                                 const VectorIntegerType& m,
                                                 BatchedDgemmComplexType** n,
                                                 const VectorIntegerType& o)
{
//	IntegerType* cptr = const_cast<IntegerType*>(&(c[0]));
//	IntegerType* dptr = const_cast<IntegerType*>(&(d[0]));
//	IntegerType* mptr = const_cast<IntegerType*>(&(m[0]));
//	IntegerType* optr = const_cast<IntegerType*>(&(o[0]));
//	setup_vbatchCmplx(a, b, cptr, dptr, e, f, g, h, i, j, k, l, mptr, n, optr);
	err("PluginSc does not support complex\n");
}

template<>
inline void setupVbatch<FpType>(IntegerType a,
                                IntegerType b,
                                const VectorIntegerType& c,
                                const VectorIntegerType& d,
                                IntegerType** e,
                                IntegerType** f,
                                IntegerType** g,
                                FpType** h,
                                IntegerType** i,
                                FpType** j,
                                IntegerType** k,
                                FpType** l,
                                const VectorIntegerType& m,
                                FpType** n,
                                const VectorIntegerType& o)
{
	IntegerType* cptr = const_cast<IntegerType*>(&(c[0]));
	IntegerType* dptr = const_cast<IntegerType*>(&(d[0]));
	IntegerType* mptr = const_cast<IntegerType*>(&(m[0]));
	IntegerType* optr = const_cast<IntegerType*>(&(o[0]));
	setup_vbatch(a, b, cptr, dptr, e, f, g, h, i, j, k, l, mptr, n, optr);
}

template<typename T>
void unsetupVbatch(IntegerType**,
                   IntegerType**,
                   IntegerType**,
                   T**,
                   IntegerType**,
                   T**,
                   IntegerType**);

template<>
inline void unsetupVbatch<FpType>(IntegerType** a,
                                  IntegerType** b,
                                  IntegerType** c,
                                  FpType** d,
                                  IntegerType** e,
                                  FpType** f,
                                  IntegerType** g)
{
	unsetup_vbatch(a, b, c, d, e, f, g);
}

template<>
inline void unsetupVbatch<BatchedDgemmComplexType>(IntegerType** a,
                                                   IntegerType** b,
                                                   IntegerType** c,
                                                   BatchedDgemmComplexType** d,
                                                   IntegerType** e,
                                                   BatchedDgemmComplexType** f,
                                                   IntegerType** g)
{
//	unsetup_vbatchCmplx(a, b, c, d, e, f, g);
	err("PluginSc does not support complex\n");
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
		VectorIntegerType ldAptr(total*total*nC);
		VectorIntegerType ldBptr(total*total*nC);

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
					aptr[outPatch + inPatch*total + ic*total*total] = a;
					bptr[outPatch + inPatch*total + ic*total*total] = b;

					ldAptr[outPatch + inPatch*total + ic*total*total] = AmatDense.rows();
					ldBptr[outPatch + inPatch*total + ic*total*total] = BmatDense.rows();
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

	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType& vout, const VectorType& vin) const
	{
		ComplexOrRealType* vinptr = const_cast<ComplexOrRealType*>(&(vin[0]));
		ComplexOrRealType* voutptr = const_cast<ComplexOrRealType*>(&(vout[0]));
		applyHtargetVbatch(initKron_.connections(), // number of connections (INPUT)
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
		                   vinptr,
		                   voutptr);
	}

private:

	void convertToVector(VectorIntegerType& v, const BasisType& b) const
	{
		SizeType total = b.partition() - 1;
		v.clear();
		v.resize(total + baseForIntegerVectors_, 0);
		for (SizeType i = 0; i < total; ++i)
			v[i + baseForIntegerVectors_] = b.partition(i + 1) - b.partition(i);
	}

	void convertOffsets(VectorIntegerType& v) const
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
	VectorIntegerType offsets_;
	VectorIntegerType pLeft_;
	VectorIntegerType pRight_;
	mutable ComplexOrRealType* Abatch_;
	mutable ComplexOrRealType* Bbatch_;
	mutable IntegerType* leftPatchStart_;
	mutable IntegerType* rightPatchStart_;
	mutable IntegerType* xyPatchStart_;
	mutable IntegerType* ldAbatch_;
	mutable IntegerType* ldBbatch_;
};
}
#endif // BATCHEDGEMM_H
