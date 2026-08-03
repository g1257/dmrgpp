#ifndef BATCHEDGEMMINCLUDE_HH
#define BATCHEDGEMMINCLUDE_HH
#include "DMRGConfig.h"
#if 0 // KOKKOS_BATCHED
#include "BatchedGemmKokkos.h"
#define BATCHED_GEMM ::Dmrg::BatchedGemmKokkos
#elif defined(PLUGIN_SC)
#include "BatchedGemmPluginSc.h"
#define BATCHED_GEMM ::Dmrg::BatchedGemmPluginSc
#else
#include "BatchedGemmCpu.h"
#define BATCHED_GEMM ::Dmrg::BatchedGemmCpu
#endif
#include <PsimagLite/PsimagLite.h>

namespace Dmrg {

class BatchedGemmInclude {

public:

	static void failIfNotSupported()
	{
#if defined(KOKKOS_BATCHED) || defined(PLUGIN_SC)
		return;
#endif
		err("BatchedGemm needs DMRG_BUILD_BATCHED_KOKKOS=ON or -DPLUGIN_SC\n");
	}

	static std::string info()
	{
#if 1 // KOKKOS_BATCHED
		return "KokkosKernels";
#elif defined(PLUGIN_SC)
		return "PLUGIN_SC";
#else
		return "";
#endif
	}
};
}
#endif // BATCHEDGEMMINCLUDE_HH
