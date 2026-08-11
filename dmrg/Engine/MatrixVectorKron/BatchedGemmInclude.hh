#ifndef BATCHEDGEMMINCLUDE_HH
#define BATCHEDGEMMINCLUDE_HH
#include "DMRGConfig.h"
#if defined(KOKKOS_BATCHED)
#include "BatchedGemmKokkos.h"
#define BATCHED_GEMM BatchedGemmKokkos
#elif defined(PLUGIN_SC)
#include "BatchedGemmPluginSc.h"
#define BATCHED_GEMM BatchedGemmPluginSc
#else
#include "BatchedGemmCpu.h"
#define BATCHED_GEMM BatchedGemmCpu
#endif
#include <PsimagLite/PsimagLite.h>

namespace Dmrg {

class BatchedGemmInclude {

public:

	static void failIfNotSupported()
	{
	}

	static std::string info()
	{
#if defined(KOKKOS_BATCHED)
		return "KOKKOS_BATCHED";
#elif defined(PLUGIN_SC)
		return "PLUGIN_SC";
#else
		return "CPU";
#endif
	}
};
}
#endif // BATCHEDGEMMINCLUDE_HH
