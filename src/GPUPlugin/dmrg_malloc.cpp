#include "dmrg_vbatch.h"
#include <cstring>
#ifdef USE_MAGMA
#include "cuda.h"
#include "cuda_runtime.h"
#endif

#include "dmrg_lapack.h"

IntegerType dmrg_is_managed(const void* ptr)
{
	const IntegerType lfalse     = (0 == 1);
	IntegerType       is_managed = lfalse;

#ifdef USE_MAGMA
	struct cudaPoIntegerTypeerAttributes attribute;
	cudaError_t ierr = cudaPoIntegerTypeerGetAttributes(&attribute, ptr);

#if defined(CUDART_VERSION) && (CUDART_VERSION >= 10000)
	is_managed = (ierr == cudaSuccess) && (attribute.type == cudaMemoryTypeManaged);
#else
	is_managed = (ierr == cudaSuccess) && attribute.isManaged;
#endif
#endif
	return (is_managed);
}

#ifdef USE_MAGMA
void dmrg_set_readonly(void* devPtr, size_t nbytes, IntegerType device)
{
	if (dmrg_is_managed(devPtr)) {

		/* cudaMemoryAdvise advice = cudaMemAdviseSetReadMostly; */
		cudaError_t istat
		    = cudaMemAdvise(devPtr, nbytes, cudaMemAdviseSetReadMostly, device);
		if (istat != cudaSuccess) {
			fprintf(stderr, "dmrg_set_readonly: %s\n", cudaGetErrorString(istat));
		};
		assert(istat == cudaSuccess);
	};
}

void dmrg_unset_readonly(void* devPtr, size_t nbytes, IntegerType device)
{
	if (dmrg_is_managed(devPtr)) {

		/* cudaMemoryAdvise advice = cudaMemAdviseUnsetReadMostly; */
		cudaError_t istat
		    = cudaMemAdvise(devPtr, nbytes, cudaMemAdviseUnsetReadMostly, device);
		if (istat != cudaSuccess) {
			fprintf(stderr, "dmrg_set_readonly: %s\n", cudaGetErrorString(istat));
		};
		assert(istat == cudaSuccess);
	};
}

#endif

void dmrg_memcpy(void* dest, const void* src, size_t count)
{
#ifdef USE_MAGMA
	cudaError_t istat = cudaMemcpy(dest, src, count, cudaMemcpyDefault);
	if (istat != cudaSuccess) {
		fprintf(stderr, "dmrg_memcpy: %s\n", cudaGetErrorString(istat));
	};
	assert(istat == cudaSuccess);
#else
	memcpy(dest, src, count);
#endif
}

template <typename T> T* dmrg_malloc(const size_t alloc_size, SizeType size)
{
	T* ptr = NULL;
#ifdef USE_MAGMA
	{
		const unsigned IntegerType flags  = cudaMemAttachGlobal;
		const size_t               nbytes = (alloc_size > 0) ? alloc_size : 1;
		cudaError_t                ierr   = cudaMallocManaged(&ptr, nbytes, flags);
		IntegerType                isok   = (ierr == cudaSuccess);
		if (!isok) {
			fprintf(stderr, "dmrg_malloc: CUDA ERROR %s\n", cudaGetErrorString(ierr));

			if (ierr == cudaErrorMemoryAllocation) {
				fprintf(stderr,
				        "dmrg_malloc:cudaErrorMemoryAllocation, alloc_size=%ld\n",
				        alloc_size);
			} else if (ierr == cudaErrorNotSupported) {
				fprintf(stderr,
				        "dmrg_malloc:cudaErrorNotSupported, alloc_size=%ld\n",
				        alloc_size);
			} else if (ierr == cudaErrorInvalidValue) {
				fprintf(stderr,
				        "dmrg_malloc:cudaErrorInvalidValue, alloc_size=%ld\n",
				        alloc_size);
			};
		};

		assert(ierr == cudaSuccess);
	}

	if (dmrg_get_ngpu() == 1) {
#ifdef USE_CUMEMADVISE
		CUdeviceptr  devPtr = (CUdeviceptr)ptr;
		size_t       count  = alloc_size;
		CUmem_advise advice = CU_MEM_ADVISE_SET_PREFERRED_LOCATION;
		CUdevice     device = 0;
		CUresult     istat  = cuMemAdvise(devPtr, count, advice, device);
		assert(istat == CUDA_SUCCESS);
#else
		IntegerType device = 0;
		size_t      count  = alloc_size;
		cudaError_t ierr
		    = cudaMemAdvise(ptr, count, cudaMemAdviseSetPreferredLocation, device);
		assert(ierr == cudaSuccess);
#endif
	}

#else
	// ptr = (void *) malloc( alloc_size );
	ptr = new T[size];
#endif
	assert(ptr != 0);
	return (ptr);
}

template <typename T>
void dmrg_lacpy(const char*       uplo,
                const IntegerType m,
                const IntegerType n,
                const T*          src,
                const IntegerType ld_src,
                T*                dest,
                const IntegerType ld_dest)
{
#ifdef USE_MAGMA
	const IntegerType is_upper = (uplo == 'U') || (uplo == 'u');
	const IntegerType is_lower = (uplo == 'L') || (uplo == 'l');
	const IntegerType is_full  = (!is_upper) && (!is_lower);

	IntegerType is_block_copy = is_full && (m == ld_src) && (m == ld_dest);
	if (is_block_copy) {
		const size_t nbytes = sizeof(T) * m * n;
		dmrg_memcpy(dest, src, nbytes);
	} else {
		const IntegerType min_mn = (m <= n) ? m : n;
		const IntegerType ncol   = is_full ? n : min_mn;
		IntegerType       jcol   = 0;

		for (jcol = 0; jcol < ncol; jcol++) {
			const IntegerType irow   = jcol;
			IntegerType       istart = is_upper ? 0 : is_lower ? irow : 0;
			IntegerType       iend   = is_upper ? irow : is_lower ? m - 1 : m - 1;
			IntegerType       count  = iend - istart + 1;
			if (count >= 1) {

				const T* psrc  = src + jcol * ld_src + istart;
				T*       pdest = dest + jcol * ld_dest + istart;

				const size_t nbytes = count * sizeof(T);
				dmrg_memcpy(pdest, psrc, nbytes);
			};
		};
	};

#else
	Xlacpy_(uplo, &m, &n, src, &ld_src, dest, &ld_dest);
#endif
}

void dmrg_prefetch_to_device(void* unified_memory_ptr, size_t nbytes, IntegerType idevice)
{
	assert(unified_memory_ptr != 0);
	if (nbytes <= 0) {
		return;
	};

#ifdef USE_MAGMA
	if (dmrg_is_managed(unified_memory_ptr)) {
		cudaError_t istat = cudaSuccess;

		IntegerType ndevice = 0;
		istat               = cudaGetDeviceCount(&ndevice);
		assert(istat == cudaSuccess);

		if (idevice > (ndevice - 1)) {
			idevice = (ndevice - 1);
		};
		if (idevice < 0) {
			idevice = 0;
		};

		istat = cudaSetDevice(idevice);
		assert(istat == cudaSuccess);

		IntegerType deviceId = 0;
		istat                = cudaGetDevice(&deviceId);
		assert(istat == cudaSuccess);

		struct cudaDeviceProp p;
		istat = cudaGetDeviceProperties(&p, deviceId);
		assert(istat == cudaSuccess);

		if (p.concurrentManagedAccess) {
			cudaStream_t stream = 0;
			const void*  devPtr = unified_memory_ptr;
			istat = cudaMemPrefetchAsync(devPtr, nbytes, deviceId, stream);
			assert(istat == cudaSuccess);

			istat = cudaDeviceSynchronize();
			assert(istat == cudaSuccess);
		}
	};
#endif
}

template <typename T> void dmrg_free(T* ptr)
{
#ifdef USE_MAGMA
	{
		cudaError_t ierr = cudaFree(ptr);
		IntegerType isok = (ierr == cudaSuccess);
		if (!isok) {
			fprintf(stderr, "dmrg_free: CUDA ERROR %s\n", cudaGetErrorString(ierr));

			if (ierr == cudaErrorInvalidDevicePoIntegerTypeer) {
				fprintf(stderr, "dmrg_free: cudaErrorInvalidDevicePointer \n");
			} else if (ierr == cudaErrorIllegalAddress) {
				fprintf(stderr, "dmrg_free: cudaErrorIllegalAddress \n");
			} else if (ierr == cudaErrorInitializationError) {
				fprintf(stderr, "dmrg_free: cudaErrorInitializationError \n");
			};
		};

		assert(ierr == cudaSuccess);
	};
#else
	//  free( ptr );
	delete ptr;
#endif
}

template MYTYPE* dmrg_malloc<MYTYPE>(const size_t alloc_size, SizeType size);
template void    dmrg_free<MYTYPE>(MYTYPE*);
template void    dmrg_lacpy<MYTYPE>(const char*       uplo,
                                 const IntegerType m,
                                 const IntegerType n,
                                 const MYTYPE*     src,
                                 const IntegerType ld_src,
                                 MYTYPE*           dest,
                                 const IntegerType ld_dest);

#if defined(USE_COMPLEX_Z)
template void    dmrg_free<double>(double*);
template double* dmrg_malloc<double>(const size_t alloc_size, SizeType size);
template void    dmrg_lacpy<double>(const char*       uplo,
                                 const IntegerType m,
                                 const IntegerType n,
                                 const double*     src,
                                 const IntegerType ld_src,
                                 double*           dest,
                                 const IntegerType ld_dest);
#else
template void                  dmrg_free<std::complex<double>>(std::complex<double>*);
template std::complex<double>* dmrg_malloc<std::complex<double>>(const size_t alloc_size,
                                                                 SizeType     size);
template void                  dmrg_lacpy<std::complex<double>>(const char*                 uplo,
                                               const IntegerType           m,
                                               const IntegerType           n,
                                               const std::complex<double>* src,
                                               const IntegerType           ld_src,
                                               std::complex<double>*       dest,
                                               const IntegerType           ld_dest);
#endif
