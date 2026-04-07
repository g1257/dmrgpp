#ifndef DMRG_VBATCH_H
#define DMRG_VBATCH_H

#include <assert.h>
#include "Vector.h"
#include "dmrg_types.h"

#ifdef USE_INTEL_MKL
#include "dmrg_mkl.h"
#endif

double dmrg_get_wtime();

typedef PsimagLite::Vector<IntegerType>::Type VectorIntegerType;
typedef typename PsimagLite::Vector<char>::Type VectorCharType;
typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

SizeType MIN(const SizeType& x,const SizeType& y);
SizeType MAX(const SizeType& x,const SizeType& y);
SizeType indx2f(const SizeType& i,const SizeType& j,const SizeType& lda);
SizeType ICEIL(const SizeType& x,const SizeType& n);
IntegerType ICEIL(const IntegerType& x,const IntegerType& n);
SizeType index3(const SizeType& ipatch,const SizeType& jpatch,const SizeType& ioperator,const SizeType& npatches);

void dmrg_init();

void dmrg_finalize();

template<typename T>
T* dmrg_malloc(const size_t alloc_size, SizeType size);

void dmrg_memcpy(void *dest, const void *src, SizeType n);

IntegerType dmrg_is_managed( const void *ptr );

#ifdef USE_MAGMA

IntegerType dmrg_get_ngpu();

IntegerType dmrg_set_max_gpus(const IntegerType max_gpus);

void dmrg_set_readonly( void *devPtr, SizeType nbytes, IntegerType device);

void dmrg_unset_readonly( void *devPtr, SizeType nbytes, IntegerType device);

void dmrg_prefetch_to_device( void *unified_memory_ptr, SizeType nbytes, IntegerType idevice );

#endif

template<typename T>	
void dmrg_lacpy(  const char* uplo,
                  const IntegerType m, 
		  const IntegerType n,
                  const T *src, 
		  const IntegerType ld_src,
                  T *dest, 
		  const IntegerType ld_dest);


template<typename T>
void  dmrg_free(T* a_ptr );

template<typename T>	
void dmrg_Xgetvector( const IntegerType n, 
                      const T *dx_src, const IntegerType incx,
                      T *hy_dst, const IntegerType incy );


template<typename T>
void dmrg_Xsetvector( const IntegerType n, 
                      const T *hx_src, const IntegerType incx,
                      T *dy_dst, const IntegerType incy );


template<typename T>	
void dmrg_Xgetmatrix( const IntegerType m, const IntegerType n, 
                      const T *dA_src, const IntegerType ldda,
                      T *hB_dst, const IntegerType ldb );


template<typename T>	
void dmrg_Xsetmatrix( const IntegerType m, const IntegerType n, 
                      const T *hA_src, const IntegerType lda,
                      T *dB_dst, const IntegerType lddb );

/*
template<typename T>	
void apply_Htarget_vbatch( 
					IntegerType noperator,
					IntegerType npatches,
					IntegerType left_patch_start_[],
					IntegerType right_patch_start_[],
					IntegerType xy_patch_start_[],
					T Abatch_[],
			IntegerType ld_Abatch,
					T Bbatch_[],
			IntegerType ld_Bbatch,
					T vin_[],
					T vout_[]);
*/

template<typename T>	
void apply_Htarget_sparse( 
        SizeType noperator,
        SizeType npatches,
        VectorSizeType left_patch_start_,
        VectorSizeType right_patch_start_,
        VectorSizeType xy_patch_start_,
        VectorSizeType nC_,
        T** gAbatch_,
        const VectorIntegerType& ld_gAbatch,
        T** gBbatch_,
        const VectorIntegerType& ld_gBbatch,
        T* X_,
        T* Y_);

template<typename T>
void dmrg_Xgemm_vbatch( char* ctransa_array,
                        char* ctransb_array,
                        IntegerType*  m_array,
                        IntegerType*  n_array,
                        IntegerType*  k_array,
                        T* alpha_array,
                        T** a_array,
                        IntegerType*    lda_array,
                        T** b_array,
                        IntegerType*    ldb_array,
                        T* beta_array,
                        T** c_array,
                        IntegerType* ldc_array,
                        SizeType group_count,
                        IntegerType* group_size);
#endif
