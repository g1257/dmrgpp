#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <complex>
#include <stdexcept>
#include <type_traits>

#include "KokkosGemm.h"
#include "KokkosType.h"

template <typename Scalar, typename IntegerForBlasType>
inline void PsimagLite::kokkos_gemm(char               transa,
                                    char               transb,
                                    IntegerForBlasType m,
                                    IntegerForBlasType n,
                                    IntegerForBlasType k,
                                    const Scalar&      alpha,
                                    const Scalar*      A,
                                    IntegerForBlasType lda,
                                    const Scalar*      B,
                                    IntegerForBlasType ldb,
                                    const Scalar&      beta,
                                    Scalar*            C,
                                    IntegerForBlasType ldc)
{
	Kokkos::Profiling::ScopedRegion scoped_region("PsimagLite::kokkos_gemm");
	int                             M      = static_cast<int>(m);
	int                             N      = static_cast<int>(n);
	int                             K      = static_cast<int>(k);
	int                             ldaVal = static_cast<int>(lda);
	int                             ldbVal = static_cast<int>(ldb);
	int                             ldcVal = static_cast<int>(ldc);

	// Normalize trans flags
	char ta = transa ? transa : 'N';
	char tb = transb ? transb : 'N';
	if (ta >= 'a' && ta <= 'z')
		ta = char(ta - 'a' + 'A');
	if (tb >= 'a' && tb <= 'z')
		tb = char(tb - 'a' + 'A');

	int req_lda = (ta == 'N') ? std::max(1, M) : std::max(1, K);
	int req_ldb = (tb == 'N') ? std::max(1, K) : std::max(1, N);
	int req_ldc = std::max(1, M);
	if (ldaVal < req_lda || ldbVal < req_ldb || ldcVal < req_ldc) {
		throw std::runtime_error("kokkos_gemm: invalid leading dimension");
	}

	// Determine Kokkos scalar type
	using KokkosScalar = KokkosType<Scalar>::type;

	Kokkos::DefaultExecutionSpace exec;
	decltype(exec)::memory_space  mem;

	// allow padded leading dimensions (ldaVal/ldbVal/ldcVal >= required)
	if (ldaVal < req_lda || ldbVal < req_ldb || ldcVal < req_ldc) {
		throw std::runtime_error("kokkos_gemm: invalid leading dimension");
	}

	using Pair = Kokkos::pair<int, int>;

	// Create host unmanaged views that reflect the actual storage (use lda/ldb as the first
	// extent) and create subviews representing the logical matrix sizes.
	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	     Aview_op(reinterpret_cast<const KokkosScalar*>(A), ldaVal, (ta == 'N' ? K : M));
	auto Aop_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Aview_op);
	auto Aop_logical = (ta == 'N') ? Kokkos::subview(Aop_device, Pair(0, M), Pair(0, K))
	                               : Kokkos::subview(Aop_device, Pair(0, K), Pair(0, M));

	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	     Bview_op(reinterpret_cast<const KokkosScalar*>(B), ldbVal, (tb == 'N' ? N : K));
	auto Bop_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Bview_op);
	auto Bop_logical = (tb == 'N') ? Kokkos::subview(Bop_device, Pair(0, K), Pair(0, N))
	                               : Kokkos::subview(Bop_device, Pair(0, N), Pair(0, K));

	// Create C view that reflects storage with possible padding (ldcVal >= M)
	Kokkos::View<KokkosScalar**, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
	     Cview(reinterpret_cast<KokkosScalar*>(C), ldcVal, N);
	auto Cop_device = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Cview);
	auto Cop_logical = Kokkos::subview(Cop_device, Pair(0, M), Pair(0, N));

	const char transA2[2] = { ta, '\0' };
	const char transB2[2] = { tb, '\0' };
	KokkosBlas::gemm(
	    exec, transA2, transB2, alpha, Aop_logical, Bop_logical, beta, Cop_logical);
	Kokkos::deep_copy(exec, Cview, Cop_device);
	exec.fence();
}

#define PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(SCALAR, INTEGER)                                        \
	template void PsimagLite::kokkos_gemm(char          transa,                                \
	                                      char          transb,                                \
	                                      INTEGER       m,                                     \
	                                      INTEGER       n,                                     \
	                                      INTEGER       k,                                     \
	                                      const SCALAR& alpha,                                 \
	                                      const SCALAR* A,                                     \
	                                      INTEGER       lda,                                   \
	                                      const SCALAR* B,                                     \
	                                      INTEGER       ldb,                                   \
	                                      const SCALAR& beta,                                  \
	                                      SCALAR*       C,                                     \
	                                      INTEGER       ldc)

#ifndef PSI_BLAS_64
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(double, int);
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(std::complex<double>, int);
#else
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(double, long int);
PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM(std::complex<double>, long int);
#endif

#undef PSIMAGLITE_INSTANTIATE_KOKKOS_GEMM
