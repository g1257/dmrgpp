#include "GPUPluginConfig.h"
#include "dmrg_vbatch.h"
#include <cstring>

#include "cuda.h"
#include "cuda_runtime.h"

IntegerType dmrg_is_managed(const void* ptr)
{
	const IntegerType lfalse     = (0 == 1);
	IntegerType       is_managed = lfalse;

	struct cudaPointerAttributes attribute;
	cudaError_t                  ierr = cudaPointerGetAttributes(&attribute, ptr);
	is_managed = (ierr == cudaSuccess) && (attribute.type == cudaMemoryTypeManaged);
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
		}
		assert(istat == cudaSuccess);
	}
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
		}
		assert(istat == cudaSuccess);
	}
}
#endif

void dmrg_memcpy(void* dest, const void* src, size_t count)
{
	cudaError_t istat = cudaMemcpy(dest, src, count, cudaMemcpyDefault);
	if (istat != cudaSuccess) {
		fprintf(stderr, "dmrg_memcpy: %s\n", cudaGetErrorString(istat));
	}
	assert(istat == cudaSuccess);
}

void dmrg_prefetch_to_device(void* unified_memory_ptr, size_t nbytes, IntegerType idevice)
{
	assert(unified_memory_ptr != 0);
	if (nbytes <= 0) {
		return;
	}
	if (dmrg_is_managed(unified_memory_ptr)) {
		IntegerType ndevice = 0;
		cudaError_t istat   = cudaGetDeviceCount(&ndevice);
		assert(istat == cudaSuccess);

		if (idevice > (ndevice - 1)) {
			idevice = (ndevice - 1);
		}
		if (idevice < 0) {
			idevice = 0;
		}

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
	}
}
