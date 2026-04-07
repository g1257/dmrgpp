#ifndef SETUP_SPARSE_BATCH_H
#define SETUP_SPARSE_BATCH_H 

#include "dmrg_types.h"

template<typename T>
void setup_sparse_batch( 
        SizeType noperator,    /* number of connections (INPUT) */
        SizeType npatches,     /* number of patches (INPUT) */
        /*
					        ----------
					        Intent(in)
					        ----------
					        */
        const VectorSizeType& left_patch_size_, /* sizes of left patches (INPUT) */
        const VectorSizeType& right_patch_size_, /* sizes of right patches (INPUT) */
        const std::vector<T*>& Amatrix_,  /* array of pointers to small Amat matrices */
        const VectorSizeType& ld_Amatirx_,   /* array of leading dimension of the small Amat matrices */
        const std::vector<T*>& Bmatrix_,  /* array of pointers to small Bmat matrices */
        const VectorSizeType& ld_Bmatrix_,  /* array of leading dimension of the small Bmat matrices */
        /*
					        -----------
					        Intent(out)
					        -----------
					        */
        VectorSizeType& pleft_patch_start, /* cumulative sum  of left_patch_size (OUTPUT) */
        VectorSizeType& pright_patch_start, /* cumulative sum of right_patch_size (OUTPUT) */
        VectorSizeType& pxy_patch_start, /* cumulative sum of product patch size (OUTPUT) */
        VectorSizeType& pnC,
        T*** pgAbatch,
        VectorIntegerType& pld_gAbatch,
        T*** pgBbatch,
        VectorIntegerType& pld_gBbatch);

template<typename T>
void unsetup_sparse_batch(
	T***,
	T***);

#endif
