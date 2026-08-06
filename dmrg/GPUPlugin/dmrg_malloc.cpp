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
	return (is_managed);
}

void dmrg_memcpy(void* dest, const void* src, size_t count) { memcpy(dest, src, count); }

template <typename T>
void dmrg_lacpy(const char*       uplo,
                const IntegerType m,
                const IntegerType n,
                const T*          src,
                const IntegerType ld_src,
                T*                dest,
                const IntegerType ld_dest)
{
	Xlacpy_(uplo, &m, &n, src, &ld_src, dest, &ld_dest);
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
