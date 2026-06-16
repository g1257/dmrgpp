#ifndef DMRG_MAGMA_H
#define DMRG_MAGMA_H 1

#include "cuda.h"
#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
#include "driver_types.h"

#include "magma_operators.h"
#include "magma_types.h"
#include "magma_v2.h"

#if defined(USE_COMPLEX_Z)

#define MAGMA_T magmaDoubleComplex

#elif defined(USE_COMPLEX_C)

#define MAGMA_T magmaFloatComplex

#elif defined(USE_FLOAT)

#define MAGMA_T float

#else

#define MAGMA_T double

#endif

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

static void magmablas_Xgemm_vbatched_max_nocheck(magma_trans_t                    transA,
                                                 magma_trans_t                    transB,
                                                 magma_IntegerType_t*             m,
                                                 magma_IntegerType_t*             n,
                                                 magma_IntegerType_t*             k,
                                                 magmaDoubleComplex               alpha,
                                                 magmaDoubleComplex const* const* dA_array,
                                                 magma_IntegerType_t*             ldda,
                                                 magmaDoubleComplex const* const* dB_array,
                                                 magma_IntegerType_t*             lddb,
                                                 magmaDoubleComplex               beta,
                                                 magmaDoubleComplex**             dC_array,
                                                 magma_IntegerType_t*             lddc,
                                                 magma_IntegerType_t              batchCount,
                                                 magma_IntegerType_t              max_m,
                                                 magma_IntegerType_t              max_n,
                                                 magma_IntegerType_t              max_k,
                                                 magma_queue_t                    queue)
{

	magmablas_zgemm_vbatched_max_nocheck(transA,
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

static void magmablas_Xgemm_vbatched_max_nocheck(magma_trans_t                   transA,
                                                 magma_trans_t                   transB,
                                                 magma_IntegerType_t*            m,
                                                 magma_IntegerType_t*            n,
                                                 magma_IntegerType_t*            k,
                                                 magmaFloatComplex               alpha,
                                                 magmaFloatComplex const* const* dA_array,
                                                 magma_IntegerType_t*            ldda,
                                                 magmaFloatComplex const* const* dB_array,
                                                 magma_IntegerType_t*            lddb,
                                                 magmaFloatComplex               beta,
                                                 magmaFloatComplex**             dC_array,
                                                 magma_IntegerType_t*            lddc,
                                                 magma_IntegerType_t             batchCount,
                                                 magma_IntegerType_t             max_m,
                                                 magma_IntegerType_t             max_n,
                                                 magma_IntegerType_t             max_k,
                                                 magma_queue_t                   queue)
{

	magmablas_cgemm_vbatched_max_nocheck(transA,
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

static void magmablas_Xgemm_vbatched_max(magma_trans_t                    transA,
                                         magma_trans_t                    transB,
                                         magma_IntegerType_t*             m,
                                         magma_IntegerType_t*             n,
                                         magma_IntegerType_t*             k,
                                         magmaDoubleComplex               alpha,
                                         magmaDoubleComplex const* const* dA_array,
                                         magma_IntegerType_t*             ldda,
                                         magmaDoubleComplex const* const* dB_array,
                                         magma_IntegerType_t*             lddb,
                                         magmaDoubleComplex               beta,
                                         magmaDoubleComplex**             dC_array,
                                         magma_IntegerType_t*             lddc,
                                         magma_IntegerType_t              batchCount,
                                         magma_IntegerType_t              max_m,
                                         magma_IntegerType_t              max_n,
                                         magma_IntegerType_t              max_k,
                                         magma_queue_t                    queue)
{

	magmablas_zgemm_vbatched_max(transA,
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

static void magmablas_Xgemm_vbatched_max(magma_trans_t                   transA,
                                         magma_trans_t                   transB,
                                         magma_IntegerType_t*            m,
                                         magma_IntegerType_t*            n,
                                         magma_IntegerType_t*            k,
                                         magmaFloatComplex               alpha,
                                         magmaFloatComplex const* const* dA_array,
                                         magma_IntegerType_t*            ldda,
                                         magmaFloatComplex const* const* dB_array,
                                         magma_IntegerType_t*            lddb,
                                         magmaFloatComplex               beta,
                                         magmaFloatComplex**             dC_array,
                                         magma_IntegerType_t*            lddc,
                                         magma_IntegerType_t             batchCount,
                                         magma_IntegerType_t             max_m,
                                         magma_IntegerType_t             max_n,
                                         magma_IntegerType_t             max_k,
                                         magma_queue_t                   queue)
{

	magmablas_cgemm_vbatched_max(transA,
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

static void magmablas_Xgemm_vbatched(magma_trans_t                   transA,
                                     magma_trans_t                   transB,
                                     magma_IntegerType_t*            m,
                                     magma_IntegerType_t*            n,
                                     magma_IntegerType_t*            k,
                                     magmaFloatComplex               alpha,
                                     magmaFloatComplex const* const* dA_array,
                                     magma_IntegerType_t*            ldda,
                                     magmaFloatComplex const* const* dB_array,
                                     magma_IntegerType_t*            lddb,
                                     magmaFloatComplex               beta,
                                     magmaFloatComplex**             dC_array,
                                     magma_IntegerType_t*            lddc,
                                     magma_IntegerType_t             batchCount,
                                     magma_queue_t                   queue)
{

	magmablas_cgemm_vbatched(transA,
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

static void magmablas_Xgemm_vbatched(magma_trans_t                    transA,
                                     magma_trans_t                    transB,
                                     magma_IntegerType_t*             m,
                                     magma_IntegerType_t*             n,
                                     magma_IntegerType_t*             k,
                                     magmaDoubleComplex               alpha,
                                     magmaDoubleComplex const* const* dA_array,
                                     magma_IntegerType_t*             ldda,
                                     magmaDoubleComplex const* const* dB_array,
                                     magma_IntegerType_t*             lddb,
                                     magmaDoubleComplex               beta,
                                     magmaDoubleComplex**             dC_array,
                                     magma_IntegerType_t*             lddc,
                                     magma_IntegerType_t              batchCount,
                                     magma_queue_t                    queue)
{

	magmablas_zgemm_vbatched(transA,
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

static void magma_Xgetmatrix(magma_IntegerType_t m,
                             magma_IntegerType_t n,
                             magmaFloatComplex*  dA_src,
                             magma_IntegerType_t ldda,
                             magmaFloatComplex*  hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_cgetmatrix(m, n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xgetmatrix(magma_IntegerType_t m,
                             magma_IntegerType_t n,
                             magmaDoubleComplex* dA_src,
                             magma_IntegerType_t ldda,
                             magmaDoubleComplex* hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_zgetmatrix(m, n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xsetmatrix(magma_IntegerType_t       m,
                             magma_IntegerType_t       n,
                             magmaDoubleComplex const* hA_src,
                             magma_IntegerType_t       lda,
                             magmaDoubleComplex*       dB_dst,
                             magma_IntegerType_t       lddb,
                             magma_queue_t             queue)
{

	magma_zsetmatrix(m, n, hA_src, lda, dB_dst, lddb, queue);
}

static void magma_Xsetmatrix(magma_IntegerType_t      m,
                             magma_IntegerType_t      n,
                             magmaFloatComplex const* hA_src,
                             magma_IntegerType_t      lda,
                             magmaFloatComplex*       dB_dst,
                             magma_IntegerType_t      lddb,
                             magma_queue_t            queue)
{

	magma_csetmatrix(m, n, hA_src, lda, dB_dst, lddb, queue);
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

static void magma_Xgetvector(magma_IntegerType_t n,
                             magmaFloatComplex*  dA_src,
                             magma_IntegerType_t ldda,
                             magmaFloatComplex*  hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_cgetvector(n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xgetvector(magma_IntegerType_t n,
                             magmaDoubleComplex* dA_src,
                             magma_IntegerType_t ldda,
                             magmaDoubleComplex* hB_dst,
                             magma_IntegerType_t ldb,
                             magma_queue_t       queue)
{

	magma_zgetvector(n, dA_src, ldda, hB_dst, ldb, queue);
}

static void magma_Xsetvector(magma_IntegerType_t       n,
                             magmaDoubleComplex const* hA_src,
                             magma_IntegerType_t       incx,
                             magmaDoubleComplex*       dB_dst,
                             magma_IntegerType_t       incy,
                             magma_queue_t             queue)
{

	magma_zsetvector(n, hA_src, incx, dB_dst, incy, queue);
}

static void magma_Xsetvector(magma_IntegerType_t      n,
                             magmaFloatComplex const* hA_src,
                             magma_IntegerType_t      incx,
                             magmaFloatComplex*       dB_dst,
                             magma_IntegerType_t      incy,
                             magma_queue_t            queue)
{

	magma_csetvector(n, hA_src, incx, dB_dst, incy, queue);
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

/*

  #if defined(USE_COMPLEX_Z)

    #define magma_Xsetmatrix magma_zsetmatrix
    #define magma_Xgetmatrix magma_zgetmatrix

    #define magma_Xsetvector magma_zsetvector
    #define magma_Xgetvector magma_zgetvector

    #define magmablas_Xgemm_vbatched magmablas_zgemm_vbatched
    #define magmablas_Xgemm_vbatched_max_nocheck magmablas_zgemm_vbatched_max_nocheck
    #define magmablas_Xgemm_vbatched_max magmablas_zgemm_vbatched_max


  #elif defined(USE_COMPLEX_C)

    #define magma_Xsetmatrix magma_csetmatrix
    #define magma_Xgetmatrix magma_cgetmatrix

    #define magma_Xsetvector magma_csetvector
    #define magma_Xgetvector magma_cgetvector

    #define magmablas_Xgemm_vbatched magmablas_cgemm_vbatched
    #define magmablas_Xgemm_vbatched_max_nocheck magmablas_cgemm_vbatched_max_nocheck
    #define magmablas_Xgemm_vbatched_max magmablas_cgemm_vbatched_max

  #elif defined(USE_FLOAT)
    #define magma_Xsetmatrix magma_ssetmatrix
    #define magma_Xgetmatrix magma_sgetmatrix

    #define magma_Xsetvector magma_ssetvector
    #define magma_Xgetvector magma_sgetvector

    #define magmablas_Xgemm_vbatched magmablas_sgemm_vbatched
    #define magmablas_Xgemm_vbatched_max_nocheck magmablas_sgemm_vbatched_max_nocheck
    #define magmablas_Xgemm_vbatched_max magmablas_sgemm_vbatched_max

  #else

    #define magma_Xsetmatrix magma_dsetmatrix
    #define magma_Xgetmatrix magma_dgetmatrix

    #define magma_Xsetvector magma_dsetvector
    #define magma_Xgetvector magma_dgetvector

    #define magmablas_Xgemm_vbatched magmablas_dgemm_vbatched
    #define magmablas_Xgemm_vbatched_max_nocheck magmablas_dgemm_vbatched_max_nocheck
    #define magmablas_Xgemm_vbatched_max magmablas_dgemm_vbatched_max

  #endif



#endif
*/

#endif
