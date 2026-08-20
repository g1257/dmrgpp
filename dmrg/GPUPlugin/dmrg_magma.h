#ifndef DMRG_MAGMA_H
#define DMRG_MAGMA_H 1

#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
#include "driver_types.h"

#include "magma_operators.h"
#include "magma_types.h"
#include "magma_v2.h"

using magma_IntegerType_t = magma_int_t;

/*
 * -------------------
 * polymorphism in C++
 * -------------------
 */

static void magmablas_Xgemm_vbatched_max_nocheck(magma_trans_t        transA,
                                                 magma_trans_t        transB,
                                                 magma_IntegerType_t* m,
                                                 magma_IntegerType_t* n,
                                                 magma_IntegerType_t* k,
                                                 double               alpha,
                                                 double const* const* dA_array,
                                                 magma_IntegerType_t* ldda,
                                                 double const* const* dB_array,
                                                 magma_IntegerType_t* lddb,
                                                 double               beta,
                                                 double**             dC_array,
                                                 magma_IntegerType_t* lddc,
                                                 magma_IntegerType_t  batchCount,
                                                 magma_IntegerType_t  max_m,
                                                 magma_IntegerType_t  max_n,
                                                 magma_IntegerType_t  max_k,
                                                 magma_queue_t        queue)
{

	magmablas_dgemm_vbatched_max_nocheck(transA,
	                                     transB,
	                                     m,
	                                     n,
	                                     k,
	                                     alpha,
	                                     dA_array,
	                                     ldda,
	                                     dB_array,
	                                     lddb,
	                                     beta,
	                                     dC_array,
	                                     lddc,
	                                     batchCount,
	                                     max_m,
	                                     max_n,
	                                     max_k,
	                                     queue);
}

static void magmablas_Xgemm_vbatched_max_nocheck(magma_trans_t        transA,
                                                 magma_trans_t        transB,
                                                 magma_IntegerType_t* m,
                                                 magma_IntegerType_t* n,
                                                 magma_IntegerType_t* k,
                                                 float                alpha,
                                                 float const* const*  dA_array,
                                                 magma_IntegerType_t* ldda,
                                                 float const* const*  dB_array,
                                                 magma_IntegerType_t* lddb,
                                                 float                beta,
                                                 float**              dC_array,
                                                 magma_IntegerType_t* lddc,
                                                 magma_IntegerType_t  batchCount,
                                                 magma_IntegerType_t  max_m,
                                                 magma_IntegerType_t  max_n,
                                                 magma_IntegerType_t  max_k,
                                                 magma_queue_t        queue)
{

	magmablas_sgemm_vbatched_max_nocheck(transA,
	                                     transB,
	                                     m,
	                                     n,
	                                     k,
	                                     alpha,
	                                     dA_array,
	                                     ldda,
	                                     dB_array,
	                                     lddb,
	                                     beta,
	                                     dC_array,
	                                     lddc,
	                                     batchCount,
	                                     max_m,
	                                     max_n,
	                                     max_k,
	                                     queue);
}

static void magmablas_Xgemm_vbatched_max_nocheck(magma_trans_t                         transA,
                                                 magma_trans_t                         transB,
                                                 magma_IntegerType_t*                  m,
                                                 magma_IntegerType_t*                  n,
                                                 magma_IntegerType_t*                  k,
                                                 std::complex<double>                  alpha,
                                                 Kokkos::complex<double> const* const* dA_array,
                                                 magma_IntegerType_t*                  ldda,
                                                 Kokkos::complex<double> const* const* dB_array,
                                                 magma_IntegerType_t*                  lddb,
                                                 std::complex<double>                  beta,
                                                 Kokkos::complex<double>**             dC_array,
                                                 magma_IntegerType_t*                  lddc,
                                                 magma_IntegerType_t                   batchCount,
                                                 magma_IntegerType_t                   max_m,
                                                 magma_IntegerType_t                   max_n,
                                                 magma_IntegerType_t                   max_k,
                                                 magma_queue_t                         queue)
{

	magmablas_zgemm_vbatched_max_nocheck(
	    transA,
	    transB,
	    m,
	    n,
	    k,
	    *reinterpret_cast<magmaDoubleComplex*>(&alpha),
	    reinterpret_cast<magmaDoubleComplex const* const*>(dA_array),
	    ldda,
	    reinterpret_cast<magmaDoubleComplex const* const*>(dB_array),
	    lddb,
	    *reinterpret_cast<magmaDoubleComplex*>(&beta),
	    reinterpret_cast<magmaDoubleComplex**>(dC_array),
	    lddc,
	    batchCount,
	    max_m,
	    max_n,
	    max_k,
	    queue);
}

static void magmablas_Xgemm_vbatched_max_nocheck(magma_trans_t                        transA,
                                                 magma_trans_t                        transB,
                                                 magma_IntegerType_t*                 m,
                                                 magma_IntegerType_t*                 n,
                                                 magma_IntegerType_t*                 k,
                                                 std::complex<float>                  alpha,
                                                 Kokkos::complex<float> const* const* dA_array,
                                                 magma_IntegerType_t*                 ldda,
                                                 Kokkos::complex<float> const* const* dB_array,
                                                 magma_IntegerType_t*                 lddb,
                                                 std::complex<float>                  beta,
                                                 Kokkos::complex<float>**             dC_array,
                                                 magma_IntegerType_t*                 lddc,
                                                 magma_IntegerType_t                  batchCount,
                                                 magma_IntegerType_t                  max_m,
                                                 magma_IntegerType_t                  max_n,
                                                 magma_IntegerType_t                  max_k,
                                                 magma_queue_t                        queue)
{

	magmablas_cgemm_vbatched_max_nocheck(
	    transA,
	    transB,
	    m,
	    n,
	    k,
	    *reinterpret_cast<magmaFloatComplex*>(&alpha),
	    reinterpret_cast<magmaFloatComplex const* const*>(dA_array),
	    ldda,
	    reinterpret_cast<magmaFloatComplex const* const*>(dB_array),
	    lddb,
	    *reinterpret_cast<magmaFloatComplex*>(&beta),
	    reinterpret_cast<magmaFloatComplex**>(dC_array),
	    lddc,
	    batchCount,
	    max_m,
	    max_n,
	    max_k,
	    queue);
}

static void magmablas_Xgemm_vbatched_max(magma_trans_t        transA,
                                         magma_trans_t        transB,
                                         magma_IntegerType_t* m,
                                         magma_IntegerType_t* n,
                                         magma_IntegerType_t* k,
                                         double               alpha,
                                         double const* const* dA_array,
                                         magma_IntegerType_t* ldda,
                                         double const* const* dB_array,
                                         magma_IntegerType_t* lddb,
                                         double               beta,
                                         double**             dC_array,
                                         magma_IntegerType_t* lddc,
                                         magma_IntegerType_t  batchCount,
                                         magma_IntegerType_t  max_m,
                                         magma_IntegerType_t  max_n,
                                         magma_IntegerType_t  max_k,
                                         magma_queue_t        queue)
{

	magmablas_dgemm_vbatched_max(transA,
	                             transB,
	                             m,
	                             n,
	                             k,
	                             alpha,
	                             dA_array,
	                             ldda,
	                             dB_array,
	                             lddb,
	                             beta,
	                             dC_array,
	                             lddc,
	                             batchCount,
	                             max_m,
	                             max_n,
	                             max_k,
	                             queue);
}

static void magmablas_Xgemm_vbatched_max(magma_trans_t        transA,
                                         magma_trans_t        transB,
                                         magma_IntegerType_t* m,
                                         magma_IntegerType_t* n,
                                         magma_IntegerType_t* k,
                                         float                alpha,
                                         float const* const*  dA_array,
                                         magma_IntegerType_t* ldda,
                                         float const* const*  dB_array,
                                         magma_IntegerType_t* lddb,
                                         float                beta,
                                         float**              dC_array,
                                         magma_IntegerType_t* lddc,
                                         magma_IntegerType_t  batchCount,
                                         magma_IntegerType_t  max_m,
                                         magma_IntegerType_t  max_n,
                                         magma_IntegerType_t  max_k,
                                         magma_queue_t        queue)
{

	magmablas_sgemm_vbatched_max(transA,
	                             transB,
	                             m,
	                             n,
	                             k,
	                             alpha,
	                             dA_array,
	                             ldda,
	                             dB_array,
	                             lddb,
	                             beta,
	                             dC_array,
	                             lddc,
	                             batchCount,
	                             max_m,
	                             max_n,
	                             max_k,
	                             queue);
}

static void magmablas_Xgemm_vbatched_max(magma_trans_t                      transA,
                                         magma_trans_t                      transB,
                                         magma_IntegerType_t*               m,
                                         magma_IntegerType_t*               n,
                                         magma_IntegerType_t*               k,
                                         std::complex<double>               alpha,
                                         std::complex<double> const* const* dA_array,
                                         magma_IntegerType_t*               ldda,
                                         std::complex<double> const* const* dB_array,
                                         magma_IntegerType_t*               lddb,
                                         std::complex<double>               beta,
                                         std::complex<double>**             dC_array,
                                         magma_IntegerType_t*               lddc,
                                         magma_IntegerType_t                batchCount,
                                         magma_IntegerType_t                max_m,
                                         magma_IntegerType_t                max_n,
                                         magma_IntegerType_t                max_k,
                                         magma_queue_t                      queue)
{

	magmablas_zgemm_vbatched_max(transA,
	                             transB,
	                             m,
	                             n,
	                             k,
	                             *reinterpret_cast<magmaDoubleComplex*>(&alpha),
	                             reinterpret_cast<magmaDoubleComplex const* const*>(dA_array),
	                             ldda,
	                             reinterpret_cast<magmaDoubleComplex const* const*>(dB_array),
	                             lddb,
	                             *reinterpret_cast<magmaDoubleComplex*>(&beta),
	                             reinterpret_cast<magmaDoubleComplex**>(dC_array),
	                             lddc,
	                             batchCount,
	                             max_m,
	                             max_n,
	                             max_k,
	                             queue);
}

static void magmablas_Xgemm_vbatched_max(magma_trans_t                     transA,
                                         magma_trans_t                     transB,
                                         magma_IntegerType_t*              m,
                                         magma_IntegerType_t*              n,
                                         magma_IntegerType_t*              k,
                                         std::complex<float>               alpha,
                                         std::complex<float> const* const* dA_array,
                                         magma_IntegerType_t*              ldda,
                                         std::complex<float> const* const* dB_array,
                                         magma_IntegerType_t*              lddb,
                                         std::complex<float>               beta,
                                         std::complex<float>**             dC_array,
                                         magma_IntegerType_t*              lddc,
                                         magma_IntegerType_t               batchCount,
                                         magma_IntegerType_t               max_m,
                                         magma_IntegerType_t               max_n,
                                         magma_IntegerType_t               max_k,
                                         magma_queue_t                     queue)
{

	magmablas_cgemm_vbatched_max(transA,
	                             transB,
	                             m,
	                             n,
	                             k,
	                             *reinterpret_cast<magmaFloatComplex*>(&alpha),
	                             reinterpret_cast<magmaFloatComplex const* const*>(dA_array),
	                             ldda,
	                             reinterpret_cast<magmaFloatComplex const* const*>(dB_array),
	                             lddb,
	                             *reinterpret_cast<magmaFloatComplex*>(&beta),
	                             reinterpret_cast<magmaFloatComplex**>(dC_array),
	                             lddc,
	                             batchCount,
	                             max_m,
	                             max_n,
	                             max_k,
	                             queue);
}

static void magmablas_Xgemm_vbatched(magma_trans_t                     transA,
                                     magma_trans_t                     transB,
                                     magma_IntegerType_t*              m,
                                     magma_IntegerType_t*              n,
                                     magma_IntegerType_t*              k,
                                     std::complex<float>               alpha,
                                     std::complex<float> const* const* dA_array,
                                     magma_IntegerType_t*              ldda,
                                     std::complex<float> const* const* dB_array,
                                     magma_IntegerType_t*              lddb,
                                     std::complex<float>               beta,
                                     std::complex<float>**             dC_array,
                                     magma_IntegerType_t*              lddc,
                                     magma_IntegerType_t               batchCount,
                                     magma_queue_t                     queue)
{

	magmablas_cgemm_vbatched(transA,
	                         transB,
	                         m,
	                         n,
	                         k,
	                         *reinterpret_cast<magmaFloatComplex*>(&alpha),
	                         reinterpret_cast<magmaFloatComplex const* const*>(dA_array),
	                         ldda,
	                         reinterpret_cast<magmaFloatComplex const* const*>(dB_array),
	                         lddb,
	                         *reinterpret_cast<magmaFloatComplex*>(&beta),
	                         reinterpret_cast<magmaFloatComplex**>(dC_array),
	                         lddc,
	                         batchCount,
	                         queue);
}

static void magmablas_Xgemm_vbatched(magma_trans_t                      transA,
                                     magma_trans_t                      transB,
                                     magma_IntegerType_t*               m,
                                     magma_IntegerType_t*               n,
                                     magma_IntegerType_t*               k,
                                     std::complex<double>               alpha,
                                     std::complex<double> const* const* dA_array,
                                     magma_IntegerType_t*               ldda,
                                     std::complex<double> const* const* dB_array,
                                     magma_IntegerType_t*               lddb,
                                     std::complex<double>               beta,
                                     std::complex<double>**             dC_array,
                                     magma_IntegerType_t*               lddc,
                                     magma_IntegerType_t                batchCount,
                                     magma_queue_t                      queue)
{
	magmablas_zgemm_vbatched(transA,
	                         transB,
	                         m,
	                         n,
	                         k,
	                         *reinterpret_cast<magmaDoubleComplex*>(&alpha),
	                         reinterpret_cast<magmaDoubleComplex const* const*>(dA_array),
	                         ldda,
	                         reinterpret_cast<magmaDoubleComplex const* const*>(dB_array),
	                         lddb,
	                         *reinterpret_cast<magmaDoubleComplex*>(&beta),
	                         reinterpret_cast<magmaDoubleComplex**>(dC_array),
	                         lddc,
	                         batchCount,
	                         queue);
}

static void magmablas_Xgemm_vbatched(magma_trans_t        transA,
                                     magma_trans_t        transB,
                                     magma_IntegerType_t* m,
                                     magma_IntegerType_t* n,
                                     magma_IntegerType_t* k,
                                     double               alpha,
                                     double const* const* dA_array,
                                     magma_IntegerType_t* ldda,
                                     double const* const* dB_array,
                                     magma_IntegerType_t* lddb,
                                     double               beta,
                                     double**             dC_array,
                                     magma_IntegerType_t* lddc,
                                     magma_IntegerType_t  batchCount,
                                     magma_queue_t        queue)
{

	magmablas_dgemm_vbatched(transA,
	                         transB,
	                         m,
	                         n,
	                         k,
	                         alpha,
	                         dA_array,
	                         ldda,
	                         dB_array,
	                         lddb,
	                         beta,
	                         dC_array,
	                         lddc,
	                         batchCount,
	                         queue);
}

static void magmablas_Xgemm_vbatched(magma_trans_t        transA,
                                     magma_trans_t        transB,
                                     magma_IntegerType_t* m,
                                     magma_IntegerType_t* n,
                                     magma_IntegerType_t* k,
                                     float                alpha,
                                     float const* const*  dA_array,
                                     magma_IntegerType_t* ldda,
                                     float const* const*  dB_array,
                                     magma_IntegerType_t* lddb,
                                     float                beta,
                                     float**              dC_array,
                                     magma_IntegerType_t* lddc,
                                     magma_IntegerType_t  batchCount,
                                     magma_queue_t        queue)
{

	magmablas_sgemm_vbatched(transA,
	                         transB,
	                         m,
	                         n,
	                         k,
	                         alpha,
	                         dA_array,
	                         ldda,
	                         dB_array,
	                         lddb,
	                         beta,
	                         dC_array,
	                         lddc,
	                         batchCount,
	                         queue);
}

static void magma_Xgetmatrix(magma_IntegerType_t m,
                             magma_IntegerType_t n,
                             double*             dA_src,
                             magma_IntegerType_t ldda,
                             double*             hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_dgetmatrix(m, n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xgetmatrix(magma_IntegerType_t m,
                             magma_IntegerType_t n,
                             float*              dA_src,
                             magma_IntegerType_t ldda,
                             float*              hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{
	magma_sgetmatrix(m, n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xgetmatrix(magma_IntegerType_t  m,
                             magma_IntegerType_t  n,
                             std::complex<float>* dA_src,
                             magma_IntegerType_t  ldda,
                             std::complex<float>* hB_dst,
                             magma_IntegerType_t  ldb,
                             magma_queue_t        queue)
{
	magma_cgetmatrix(m,
	                 n,
	                 reinterpret_cast<magmaFloatComplex*>(dA_src),
	                 ldda,
	                 reinterpret_cast<magmaFloatComplex*>(hB_dst),
	                 ldb,
	                 queue);
}

static void magma_Xgetmatrix(magma_IntegerType_t   m,
                             magma_IntegerType_t   n,
                             std::complex<double>* dA_src,
                             magma_IntegerType_t   ldda,
                             std::complex<double>* hB_dst,
                             magma_IntegerType_t   ldb,
                             magma_queue_t         queue)
{
	magma_zgetmatrix(m,
	                 n,
	                 reinterpret_cast<magmaDoubleComplex*>(dA_src),
	                 ldda,
	                 reinterpret_cast<magmaDoubleComplex*>(hB_dst),
	                 ldb,
	                 queue);
}

static void magma_Xsetmatrix(magma_IntegerType_t         m,
                             magma_IntegerType_t         n,
                             std::complex<double> const* hA_src,
                             magma_IntegerType_t         lda,
                             std::complex<double>*       dB_dst,
                             magma_IntegerType_t         lddb,
                             magma_queue_t               queue)
{
	magma_zsetmatrix(m,
	                 n,
	                 reinterpret_cast<const magmaDoubleComplex*>(hA_src),
	                 lda,
	                 reinterpret_cast<magmaDoubleComplex*>(dB_dst),
	                 lddb,
	                 queue);
}

static void magma_Xsetmatrix(magma_IntegerType_t        m,
                             magma_IntegerType_t        n,
                             std::complex<float> const* hA_src,
                             magma_IntegerType_t        lda,
                             std::complex<float>*       dB_dst,
                             magma_IntegerType_t        lddb,
                             magma_queue_t              queue)
{
	magma_csetmatrix(m,
	                 n,
	                 reinterpret_cast<const magmaFloatComplex*>(hA_src),
	                 lda,
	                 reinterpret_cast<magmaFloatComplex*>(dB_dst),
	                 lddb,
	                 queue);
}

static void magma_Xsetmatrix(magma_IntegerType_t m,
                             magma_IntegerType_t n,
                             float const*        hA_src,
                             magma_IntegerType_t lda,
                             float*              dB_dst,
                             magma_IntegerType_t lddb,
                             magma_queue_t       queue)
{

	magma_ssetmatrix(m, n, hA_src, lda, dB_dst, lddb, queue);
}

static void magma_Xsetmatrix(magma_IntegerType_t m,
                             magma_IntegerType_t n,
                             double const*       hA_src,
                             magma_IntegerType_t lda,
                             double*             dB_dst,
                             magma_IntegerType_t lddb,
                             magma_queue_t       queue)
{

	magma_dsetmatrix(m, n, hA_src, lda, dB_dst, lddb, queue);
}

static void magma_Xgetvector(magma_IntegerType_t n,
                             double*             dA_src,
                             magma_IntegerType_t ldda,
                             double*             hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_dgetvector(n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xgetvector(magma_IntegerType_t n,
                             float*              dA_src,
                             magma_IntegerType_t ldda,
                             float*              hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_sgetvector(n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xgetvector(magma_IntegerType_t  n,
                             std::complex<float>* dA_src,
                             magma_IntegerType_t  ldda,
                             std::complex<float>* hB_dst,
                             magma_IntegerType_t  ldb,
                             magma_queue_t        queue)
{

	magma_cgetvector(n,
	                 reinterpret_cast<magmaFloatComplex*>(dA_src),
	                 ldda,
	                 reinterpret_cast<magmaFloatComplex*>(hB_dst),
	                 ldb,
	                 queue);
}

static void magma_Xgetvector(magma_IntegerType_t   n,
                             std::complex<double>* dA_src,
                             magma_IntegerType_t   ldda,
                             std::complex<double>* hB_dst,
                             magma_IntegerType_t   ldb,
                             magma_queue_t         queue)
{

	magma_zgetvector(n,
	                 reinterpret_cast<magmaDoubleComplex*>(dA_src),
	                 ldda,
	                 reinterpret_cast<magmaDoubleComplex*>(hB_dst),
	                 ldb,
	                 queue);
}

static void magma_Xsetvector(magma_IntegerType_t         n,
                             std::complex<double> const* hA_src,
                             magma_IntegerType_t         incx,
                             std::complex<double>*       dB_dst,
                             magma_IntegerType_t         incy,
                             magma_queue_t               queue)
{

	magma_zsetvector(n,
	                 reinterpret_cast<const magmaDoubleComplex*>(hA_src),
	                 incx,
	                 reinterpret_cast<magmaDoubleComplex*>(dB_dst),
	                 incy,
	                 queue);
}

static void magma_Xsetvector(magma_IntegerType_t        n,
                             std::complex<float> const* hA_src,
                             magma_IntegerType_t        incx,
                             std::complex<float>*       dB_dst,
                             magma_IntegerType_t        incy,
                             magma_queue_t              queue)
{

	magma_csetvector(n,
	                 reinterpret_cast<const magmaFloatComplex*>(hA_src),
	                 incx,
	                 reinterpret_cast<magmaFloatComplex*>(dB_dst),
	                 incy,
	                 queue);
}

static void magma_Xsetvector(magma_IntegerType_t n,
                             float const*        hA_src,
                             magma_IntegerType_t incx,
                             float*              dB_dst,
                             magma_IntegerType_t incy,
                             magma_queue_t       queue)
{

	magma_ssetvector(n, hA_src, incx, dB_dst, incy, queue);
}

static void magma_Xsetvector(magma_IntegerType_t n,
                             double const*       hA_src,
                             magma_IntegerType_t incx,
                             double*             dB_dst,
                             magma_IntegerType_t incy,
                             magma_queue_t       queue)
{

	magma_dsetvector(n, hA_src, incx, dB_dst, incy, queue);
}

#endif
