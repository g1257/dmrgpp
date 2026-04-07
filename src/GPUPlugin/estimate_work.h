#ifndef ESTIMATE_WORK_H
#define ESTIMATE_WORK_H

#include "Vector.h"

typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

void estimate_work(SizeType              npatches,
                   const VectorSizeType& left_patch_size_,
                   const VectorSizeType& right_patch_size_,
                   VectorSizeType&       nC_,
                   double*               ptotal_gflops,
                   double*               pgmemA,
                   double*               pgmemB,
                   double*               pgmemBX,
                   double*               pgmemXY);

template <typename T>
void get_total_memory(SizeType               noperator,
                      SizeType               npatches,
                      const std::vector<T*>& Amatrix_,
                      const VectorSizeType&  ld_Amatrix_,
                      const std::vector<T*>& Bmatrix_,
                      const VectorSizeType&  ld_Bmatrix_,
                      const VectorSizeType&  left_patch_size_,
                      const VectorSizeType&  right_patch_size_,
                      SizeType&              ptotal_memory_in_nbytes);

#endif
