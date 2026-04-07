#ifndef SETUP_VBATCH_H
#define SETUP_VBATCH_H

#include "Vector.h"
#include "dmrg_types.h"

typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

#if defined(USE_COMPLEX_Z)
double ABS(std::complex<double> x);
#else
double ABS(double x);
#endif

template <typename T>
void setup_vbatch(
    SizeType               noperator, /* number of connections (INPUT) */
    SizeType               npatches, /* number of patches (INPUT) */
    const VectorSizeType&  left_patch_size_, /* sizes of left patches (INPUT) */
    const VectorSizeType&  right_patch_size_, /* sizes of right patches (INPUT) */
    VectorSizeType&        left_patch_start_, /* cumulative sum  of left_patch_size (OUTPUT) */
    VectorSizeType&        right_patch_start_, /* cumulative sum of right_patch_size (OUTPUT) */
    VectorSizeType&        xy_patch_start_, /* cumulative sum of product patch size (OUTPUT) */
    std::vector<T>&        pAbatch, /* Abatch matrix (OUTPUT) */
    SizeType&              pld_Abatch, /* leading dimension of Abatch matrix (OUTPUT) */
    std::vector<T>&        pBbatch, /* Bbatch matrix (OUTPUT) */
    SizeType&              pld_Bbatch, /*  leading dimension of Bbatch matrix (OUTPUT) */
    const std::vector<T*>& Amatrix_, /* array of pointers to small Amat matrices */
    const VectorSizeType&  ld_Amatrix_, /* array of leading dimension of the small Amat matrices */
    const std::vector<T*>& Bmatrix_, /* array of pointers to small Bmat matrices */
    const VectorSizeType&  ld_Bmatrix_ /* array of leading dimension of the small Bmat matrices */
);

#endif
