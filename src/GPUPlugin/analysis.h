#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "Vector.h"
#include <math.h>

typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

SizeType indx2f2(const SizeType& i, const SizeType& j, const SizeType& lda);
SizeType MOD(const SizeType& x, const SizeType& y);

SizeType gen_patches_comb(SizeType        left_size,
                          SizeType        right_size,
                          SizeType        target_up,
                          SizeType        target_down,
                          SizeType        keep_left_states,
                          SizeType        keep_right_states,
                          VectorSizeType& left_patch_size_,
                          VectorSizeType& right_patch_size_,
                          VectorSizeType& left_patch_up_,
                          VectorSizeType& left_patch_down_,
                          VectorSizeType& right_patch_up_,
                          VectorSizeType& right_patch_down_);

void cal_kron_flops(SizeType nrowA,
                    SizeType nrowB,
                    SizeType ncolA,
                    SizeType ncolB,
                    double*  pflops_total,
                    double*  pflops_method1,
                    double*  pflops_method2);

#endif
