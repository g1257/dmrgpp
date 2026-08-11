#ifndef SETUP_MATRIX_H
#define SETUP_MATRIX_H

#include <PsimagLite/Vector.h>

#include "dmrg_types.h"

using VectorSizeType = PsimagLite::Vector<SizeType>::Type;

template <typename T>
void setup_matrix(SizeType              noperator,
                  SizeType              npatches,
                  std::vector<T>&       Abatch,
                  const VectorSizeType& left_patch_size_,
                  std::vector<T*>&      Amatrix_,
                  VectorSizeType&       ld_Amatrix_);

#endif
