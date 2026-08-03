#ifndef PSIMAG_KOKKOS_GEMM_H
#define PSIMAG_KOKKOS_GEMM_H

namespace PsimagLite {
template <typename Scalar, typename IntegerForBlasType>
void kokkos_gemm(char               transa,
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
                 IntegerForBlasType ldc);
} // namespace PsimagLite

#endif // PSIMAG_KOKKOS_GEMM_H
