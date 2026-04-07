#ifndef SETUP_NC_H
#define SETUP_NC_H

#include "dmrg_types.h"

template <typename T>
void setup_nC(SizeType noperator,
              SizeType npatches,
              /*
               * ---------
               * Intent(in)
               * ---------
               */
              const std::vector<T*>& Amatrix_,
              const VectorSizeType&  ld_Amatrix_,
              const std::vector<T*>& Bmatrix_,
              const VectorSizeType&  ld_Bmatrix_,
              const VectorSizeType&  left_patch_size_, /* sizes of left patches (INPUT) */
              const VectorSizeType&  right_patch_size_, /* sizes of right patches (INPUT) */
              /*
               * --------------------------------
               * Intent(out) but assume storage already allocated
               * --------------------------------
               */
              VectorSizeType& nC_,
              VectorSizeType& gnnz_A_,
              VectorSizeType& gnnz_B_);

#endif
