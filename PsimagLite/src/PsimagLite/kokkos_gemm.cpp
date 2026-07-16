#include <KokkosBlas3_gemm.hpp>
#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <complex>
#include <stdexcept>
#include <type_traits>

#include <PsimagLite/kokkos_gemm.h>

namespace {

// The scalar types that are floating point types and their corresponding std::complex types.
// We need to map std::complex to Kokkos::complex while keeping all the other types which is what
// KokkosType does.
template <typename T> struct KokkosType {
	using type = T;
};

template <typename T>
        requires(!std::is_floating_point_v<T>)
struct KokkosType<T> {
	using type = Kokkos::complex<typename T::value_type>;
};

}

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

	KOKKOS_ASSERT(ldaVal == (ta == 'N' ? M : K));
	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	    Aview_op(
	        reinterpret_cast<const KokkosScalar*>(A), ta == 'N' ? M : K, ta == 'N' ? K : M);
	auto Aview_op_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Aview_op);

	KOKKOS_ASSERT(ldbVal == (tb == 'N' ? K : N));
	Kokkos::View<const KokkosScalar**,
	             Kokkos::LayoutLeft,
	             Kokkos::HostSpace,
	             Kokkos::MemoryUnmanaged>
	    Bview_op(
	        reinterpret_cast<const KokkosScalar*>(B), tb == 'N' ? K : N, tb == 'N' ? N : K);
	auto Bview_op_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Bview_op);

	Kokkos::View<KokkosScalar**, Kokkos::LayoutLeft, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>
	     Cview(reinterpret_cast<KokkosScalar*>(C), M, N);
	auto Cview_device
	    = Kokkos::create_mirror_view_and_copy(Kokkos::view_alloc(exec, mem), Cview);

	const char transA[2] = { ta, '\0' };
	const char transB[2] = { tb, '\0' };
	KokkosBlas::gemm(
	    exec, transA, transB, alpha, Aview_op_device, Bview_op_device, beta, Cview_device);
	Kokkos::deep_copy(exec, Cview, Cview_device);
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
