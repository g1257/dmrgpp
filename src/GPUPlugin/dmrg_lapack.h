#ifndef DMRG_LAPACK_H
#define DMRG_LAPACK_H

#include "dmrg_types.h"
#include <complex>

/*
 -----------------------------
 polymorphism supported in C++
 -----------------------------
 */

extern "C"
{
	extern void zlacpy_(const char*                 uplo,
	                    const IntegerType*          m,
	                    const IntegerType*          n,
	                    const std::complex<double>* A,
	                    const IntegerType*          lda,
	                    std::complex<double>*       B,
	                    const IntegerType*          ldb);
	extern void dlacpy_(const char*        uplo,
	                    const IntegerType* m,
	                    const IntegerType* n,
	                    const double*      A,
	                    const IntegerType* lda,
	                    double*            B,
	                    const IntegerType* ldb);
	extern void zgemm_(const char*                 transa,
	                   const char*                 transb,
	                   const IntegerType*          m,
	                   const IntegerType*          n,
	                   const IntegerType*          k,
	                   const std::complex<double>* alpha,
	                   const std::complex<double>* A,
	                   const IntegerType*          lda,
	                   const std::complex<double>* B,
	                   const IntegerType*          ldb,
	                   const std::complex<double>* beta,
	                   std::complex<double>*       C,
	                   const IntegerType*          ldc);
	extern void dgemm_(const char*        transa,
	                   const char*        transb,
	                   const IntegerType* m,
	                   const IntegerType* n,
	                   const IntegerType* k,
	                   const double*      alpha,
	                   const double*      A,
	                   const IntegerType* lda,
	                   const double*      B,
	                   const IntegerType* ldb,
	                   const double*      beta,
	                   double*            C,
	                   const IntegerType* ldc);
	extern void zcopy_(const IntegerType*          n,
	                   const std::complex<double>* x,
	                   const IntegerType*          incx,
	                   std::complex<double>*       y,
	                   const IntegerType*          incy);
	extern void dcopy_(const IntegerType* n,
	                   const double*      x,
	                   const IntegerType* incx,
	                   double*            y,
	                   const IntegerType* incy);
}

void Xlacpy_(const char*                 uplo,
             const IntegerType*          m,
             const IntegerType*          n,
             const std::complex<double>* A,
             const IntegerType*          lda,
             std::complex<double>*       B,
             const IntegerType*          ldb);

void Xgemm_(const char*                 transA,
            const char*                 transB,
            const IntegerType*          m,
            const IntegerType*          n,
            const IntegerType*          k,
            const std::complex<double>* alpha,
            const std::complex<double>* A,
            const IntegerType*          lda,
            const std::complex<double>* B,
            const IntegerType*          ldb,
            const std::complex<double>* beta,
            std::complex<double>*       C,
            const IntegerType*          ldc);

void Xcopy_(const IntegerType*          n,
            const std::complex<double>* X,
            const IntegerType*          incX,
            std::complex<double>*       Y,
            const IntegerType*          incY);

void Xlacpy_(const char*        uplo,
             const IntegerType* m,
             const IntegerType* n,
             const double*      A,
             const IntegerType* lda,
             double*            B,
             const IntegerType* ldb);

void Xgemm_(const char*        transA,
            const char*        transB,
            const IntegerType* m,
            const IntegerType* n,
            const IntegerType* k,
            const double*      alpha,
            const double*      A,
            const IntegerType* lda,
            const double*      B,
            const IntegerType* ldb,
            const double*      beta,
            double*            C,
            const IntegerType* ldc);

void Xcopy_(const IntegerType* n,
            const double*      X,
            const IntegerType* incX,
            double*            Y,
            const IntegerType* incY);

#endif
