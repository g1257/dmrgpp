#include "GPUPluginConfig.h"
#include "dmrg_vbatch.h"
#include <cstring>

#include "hip/hip_runtime.h"

IntegerType dmrg_is_managed(const void* ptr)
{
	const IntegerType lfalse     = (0 == 1);
	IntegerType       is_managed = lfalse;

	struct hipPointerAttribute_t attribute;
	hipError_t                   ierr = hipPointerGetAttributes(&attribute, ptr);

	is_managed = (ierr == hipSuccess) && (attribute.type == hipMemoryTypeManaged);

	return (is_managed);
}

#ifdef USE_MAGMA
void dmrg_set_readonly(void* devPtr, size_t nbytes, IntegerType device)
{
	if (dmrg_is_managed(devPtr)) {

		/* cudaMemoryAdvise advice = cudaMemAdviseSetReadMostly; */
		hipMemLocation gpuLocation;
		gpuLocation.type = hipMemLocationTypeDevice; // or cudaMemLocationTypeHostNuma
		gpuLocation.id   = device;
		hipError_t istat
		    = hipMemAdvise_v2(devPtr, nbytes, hipMemAdviseSetReadMostly, gpuLocation);
		if (istat != hipSuccess) {
			fprintf(stderr, "dmrg_set_readonly: %s\n", hipGetErrorString(istat));
		}
		assert(istat == hipSuccess);
	}
}

void dmrg_unset_readonly(void* devPtr, size_t nbytes, IntegerType device)
{
	if (dmrg_is_managed(devPtr)) {

		/* cudaMemoryAdvise advice = cudaMemAdviseUnsetReadMostly; */
		hipMemLocation gpuLocation;
		gpuLocation.type = hipMemLocationTypeDevice; // or cudaMemLocationTypeHostNuma
		gpuLocation.id   = device;
		hipError_t istat
		    = hipMemAdvise_v2(devPtr, nbytes, hipMemAdviseUnsetReadMostly, gpuLocation);
		if (istat != hipSuccess) {
			fprintf(stderr, "dmrg_set_readonly: %s\n", hipGetErrorString(istat));
		}
		assert(istat == hipSuccess);
	}
}
#endif

void dmrg_memcpy(void* dest, const void* src, size_t count)
{
	hipError_t istat = hipMemcpy(dest, src, count, hipMemcpyDefault);
	if (istat != hipSuccess) {
		fprintf(stderr, "dmrg_memcpy: %s\n", hipGetErrorString(istat));
	}
	assert(istat == hipSuccess);
}

void dmrg_prefetch_to_device(void* unified_memory_ptr, size_t nbytes, IntegerType idevice)
{
	assert(unified_memory_ptr != 0);
	if (nbytes <= 0) {
		return;
	}
	if (dmrg_is_managed(unified_memory_ptr)) {
		IntegerType ndevice = 0;
		hipError_t  istat   = hipGetDeviceCount(&ndevice);
		assert(istat == hipSuccess);

		if (idevice > (ndevice - 1)) {
			idevice = (ndevice - 1);
		}
		if (idevice < 0) {
			idevice = 0;
		}

		istat = hipSetDevice(idevice);
		assert(istat == hipSuccess);

		IntegerType deviceId = 0;
		istat                = hipGetDevice(&deviceId);
		assert(istat == hipSuccess);

		hipMemLocation gpuLocation;
		gpuLocation.type = hipMemLocationTypeDevice; // or cudaMemLocationTypeHostNuma
		gpuLocation.id   = deviceId;

		struct hipDeviceProp_t p;
		istat = hipGetDeviceProperties(&p, deviceId);
		assert(istat == hipSuccess);

		if (p.concurrentManagedAccess) {
			hipStream_t stream;
			hipError_t  result = hipStreamCreate(&stream);
			if (result != hipSuccess) {
				// Handle error
				fprintf(stderr,
				        "hipStreamCreate failed: %s\n",
				        hipGetErrorString(result));
			}
			const void* devPtr = unified_memory_ptr;
			istat = hipMemPrefetchAsync_v2(devPtr, nbytes, gpuLocation, 0, stream);
			assert(istat == hipSuccess);

			istat = hipDeviceSynchronize();
			assert(istat == hipSuccess);
		}
#ifdef NDEBUG
		(void)istat;
#endif
	}
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
