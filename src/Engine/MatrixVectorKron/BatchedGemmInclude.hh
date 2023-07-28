#ifndef BATCHEDGEMMINCLUDE_HH
#define BATCHEDGEMMINCLUDE_HH
#ifdef PLUGIN_SC
#include "BatchedGemmPluginSc.h"
#define BATCHED_GEMM BatchedGemmPluginSc
#else
#include "BatchedGemmCpu.h"
#define BATCHED_GEMM BatchedGemmCpu
#endif

namespace Dmrg
{

class BatchedGemmInclude
{

public:

	static void failIfNotSupported()
	{
#ifdef PLUGIN_SC
		return;
#endif
		err("BatchedGemm needs -DPLUGIN_SC in Config.make\n");
	}

	static std::string info()
	{
#ifdef PLUGIN_SC
		return "PLUGIN_SC";
#else
		return "";
#endif
	}
};
}
#endif // BATCHEDGEMMINCLUDE_HH
