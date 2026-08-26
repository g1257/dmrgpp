#include "GPUPluginConfig.h"
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
	struct cudaPointerAttributes attribute;
	cudaError_t                  ierr = cudaPointerGetAttributes(&attribute, ptr);

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
		cudaMemLocation gpuLocation;
		gpuLocation.type = cudaMemLocationTypeDevice; // or cudaMemLocationTypeHostNuma
		gpuLocation.id   = device;
		cudaError_t istat
		    = cudaMemAdvise(devPtr, nbytes, cudaMemAdviseSetReadMostly, gpuLocation);
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
		cudaMemLocation gpuLocation;
		gpuLocation.type = cudaMemLocationTypeDevice; // or cudaMemLocationTypeHostNuma
		gpuLocation.id   = device;
		cudaError_t istat
		    = cudaMemAdvise(devPtr, nbytes, cudaMemAdviseUnsetReadMostly, gpuLocation);
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
	const IntegerType is_upper = (*uplo == 'U') || (*uplo == 'u');
	const IntegerType is_lower = (*uplo == 'L') || (*uplo == 'l');
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
		IntegerType ndevice = 0;
		cudaError_t istat   = cudaGetDeviceCount(&ndevice);
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

		cudaMemLocation gpuLocation;
		gpuLocation.type = cudaMemLocationTypeDevice; // or cudaMemLocationTypeHostNuma
		gpuLocation.id   = deviceId;

		struct cudaDeviceProp p;
		istat = cudaGetDeviceProperties(&p, deviceId);
		assert(istat == cudaSuccess);

		if (p.concurrentManagedAccess) {
			cudaStream_t stream;
			cudaError_t  result = cudaStreamCreate(&stream);
			if (result != cudaSuccess) {
				// Handle error
				fprintf(stderr,
				        "cudaStreamCreate failed: %s\n",
				        cudaGetErrorString(result));
			}
			const void* devPtr = unified_memory_ptr;
			istat = cudaMemPrefetchAsync(devPtr, nbytes, gpuLocation, 0, stream);
			assert(istat == cudaSuccess);

			istat = cudaDeviceSynchronize();
			assert(istat == cudaSuccess);
		}
#ifdef NDEBUG
		(void)istat;
#endif
	};
#endif
}

template void dmrg_lacpy<double>(const char*       uplo,
                                 const IntegerType m,
                                 const IntegerType n,
                                 const double*     src,
                                 const IntegerType ld_src,
                                 double*           dest,
                                 const IntegerType ld_dest);

template void dmrg_lacpy<std::complex<double>>(const char*                 uplo,
                                               const IntegerType           m,
                                               const IntegerType           n,
                                               const std::complex<double>* src,
                                               const IntegerType           ld_src,
                                               std::complex<double>*       dest,
                                               const IntegerType           ld_dest);
