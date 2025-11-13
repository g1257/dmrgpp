//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for BLAS                                                     *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// ****************************************************************************

#ifndef PSIMAG_LAPACK
#define PSIMAG_LAPACK

#include <complex>

/** \file LAPACK
 *  \author Thomas C. Schulthess,
 MSS
 */

namespace psimag {

/** \brief Namespace encapsulating psimag wrappers of LAPACK functions.
 */
namespace LAPACK {

#ifndef PSI_LAPACK_64
	typedef int IntegerForLapackType;
#else
	typedef long int IntegerForLapackType;
#endif
	// ============================================================================
	extern "C" void sgesv_(IntegerForLapackType*, IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void dgesv_(IntegerForLapackType*, IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void cgesv_(IntegerForLapackType*, IntegerForLapackType*, std::complex<float>*, IntegerForLapackType*, IntegerForLapackType*, std::complex<float>*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void zgesv_(IntegerForLapackType*, IntegerForLapackType*, std::complex<double>*, IntegerForLapackType*, IntegerForLapackType*, std::complex<double>*, IntegerForLapackType*, IntegerForLapackType*);

	// MSS
	extern "C" IntegerForLapackType
	dgetrf_(IntegerForLapackType*, IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType
	sgetrf_(IntegerForLapackType*, IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType
	zgetrf_(IntegerForLapackType*, IntegerForLapackType*, std::complex<double>*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType
	cgetrf_(IntegerForLapackType*, IntegerForLapackType*, std::complex<float>*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType
	zgetri_(IntegerForLapackType*, std::complex<double>*, IntegerForLapackType*, IntegerForLapackType*, std::complex<double>*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType
	cgetri_(IntegerForLapackType*, std::complex<float>*, IntegerForLapackType*, IntegerForLapackType*, std::complex<float>*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType dgetri_(IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType sgetri_(IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" IntegerForLapackType
	dgesdd_(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, double* a, IntegerForLapackType* lda, double* s, double* u, IntegerForLapackType* ldu, double* vt, IntegerForLapackType* ldvt, double* work, IntegerForLapackType* lwork, IntegerForLapackType* iwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	sgesdd_(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, float* a, IntegerForLapackType* lda, float* s, float* u, IntegerForLapackType* ldu, float* vt, IntegerForLapackType* ldvt, float* work, IntegerForLapackType* lwork, IntegerForLapackType* iwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	zgesdd_(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<double>* a, IntegerForLapackType* lda, double* s, std::complex<double>* u, IntegerForLapackType* ldu, std::complex<double>* vt, IntegerForLapackType* ldvt, std::complex<double>* work, IntegerForLapackType* lwork, double* rwork, IntegerForLapackType* iwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	cgesdd_(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<float>* a, IntegerForLapackType* lda, float* s, std::complex<float>* u, IntegerForLapackType* ldu, std::complex<float>* vt, IntegerForLapackType* ldvt, std::complex<float>* work, IntegerForLapackType* lwork, float* rwork, IntegerForLapackType* iwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	dgesvd_(char* jobz, char*, IntegerForLapackType* m, IntegerForLapackType* n, double* a, IntegerForLapackType* lda, double* s, double* u, IntegerForLapackType* ldu, double* vt, IntegerForLapackType* ldvt, double* work, IntegerForLapackType* lwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	sgesvd_(char* jobz, char*, IntegerForLapackType* m, IntegerForLapackType* n, float* a, IntegerForLapackType* lda, float* s, float* u, IntegerForLapackType* ldu, float* vt, IntegerForLapackType* ldvt, float* work, IntegerForLapackType* lwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	zgesvd_(char* jobz, char*, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<double>* a, IntegerForLapackType* lda, double* s, std::complex<double>* u, IntegerForLapackType* ldu, std::complex<double>* vt, IntegerForLapackType* ldvt, std::complex<double>* work, IntegerForLapackType* lwork, double* rwork, IntegerForLapackType* info);

	extern "C" IntegerForLapackType
	cgesvd_(char* jobz, char*, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<float>* a, IntegerForLapackType* lda, float* s, std::complex<float>* u, IntegerForLapackType* ldu, std::complex<float>* vt, IntegerForLapackType* ldvt, std::complex<float>* work, IntegerForLapackType* lwork, float* rwork, IntegerForLapackType* info);

	extern "C" void ilaver_(IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void dstedc_(char*, IntegerForLapackType*, double*, double*, double*, IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void sstedc_(char*, IntegerForLapackType*, float*, float*, float*, IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void zstedc_(char*, IntegerForLapackType*, double*, double*, std::complex<double>*, IntegerForLapackType*, std::complex<double>*, IntegerForLapackType*, double*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void cstedc_(char*, IntegerForLapackType*, float*, float*, std::complex<float>*, IntegerForLapackType*, std::complex<float>*, IntegerForLapackType*, float*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*, IntegerForLapackType*);

	extern "C" void dsterf_(IntegerForLapackType*, double*, double*, IntegerForLapackType*);

	extern "C" void ssterf_(IntegerForLapackType*, float*, float*, IntegerForLapackType*);

	extern "C" void dsteqr_(char*, IntegerForLapackType*, double*, double*, double*, IntegerForLapackType*, double*, IntegerForLapackType*);

	extern "C" void ssteqr_(char*, IntegerForLapackType*, float*, float*, float*, IntegerForLapackType*, float*, IntegerForLapackType*);

	// ============================================================================

	inline void STERF(IntegerForLapackType* n, double* d, double* e, IntegerForLapackType* info)
	{
		dsterf_(n, d, e, info);
	}

	inline void STERF(IntegerForLapackType* n, float* d, float* e, IntegerForLapackType* info)
	{
		ssterf_(n, d, e, info);
	}

	inline void STEQR(char jobz, IntegerForLapackType n, double* d, double* e, double* z, IntegerForLapackType ldz, double* work, IntegerForLapackType* info)
	{
		dsteqr_(&jobz, &n, d, e, z, &ldz, work, info);
	}

	inline void STEQR(char jobz, IntegerForLapackType n, float* d, float* e, float* z, IntegerForLapackType ldz, float* work, IntegerForLapackType* info)
	{
		ssteqr_(&jobz, &n, d, e, z, &ldz, work, info);
	}

	inline void GESV(IntegerForLapackType ma, IntegerForLapackType mb, float* a, IntegerForLapackType lda, IntegerForLapackType* pivot, float* b, IntegerForLapackType ldb, int& info)
	{
		sgesv_(&ma, &mb, a, &lda, pivot, b, &ldb, &info);
	}

	inline void GESV(IntegerForLapackType ma, IntegerForLapackType mb, double* a, IntegerForLapackType lda, IntegerForLapackType* pivot, double* b, IntegerForLapackType ldb, int& info)
	{
		dgesv_(&ma, &mb, a, &lda, pivot, b, &ldb, &info);
	}

	inline void GESV(IntegerForLapackType ma, IntegerForLapackType mb, std::complex<float>* a, IntegerForLapackType lda, IntegerForLapackType* pivot, std::complex<float>* b, IntegerForLapackType ldb, int& info)
	{
		cgesv_(&ma, &mb, a, &lda, pivot, b, &ldb, &info);
	}

	inline void GESV(IntegerForLapackType ma, IntegerForLapackType mb, std::complex<double>* a, IntegerForLapackType lda, IntegerForLapackType* pivot, std::complex<double>* b, IntegerForLapackType ldb, int& info)
	{
		zgesv_(&ma, &mb, a, &lda, pivot, b, &ldb, &info);
	}

	inline void GETRF(IntegerForLapackType ma, IntegerForLapackType na, double* a, IntegerForLapackType lda, IntegerForLapackType* pivot, int& info)
	{
		dgetrf_(&ma, &na, a, &lda, pivot, &info);
	}

	inline void GETRF(IntegerForLapackType ma, IntegerForLapackType na, std::complex<double>* a, IntegerForLapackType lda, IntegerForLapackType* pivot, int& info)
	{
		zgetrf_(&ma, &na, a, &lda, pivot, &info);
	}

	inline void GETRF(IntegerForLapackType ma, IntegerForLapackType na, float* a, IntegerForLapackType lda, IntegerForLapackType* pivot, int& info)
	{
		sgetrf_(&ma, &na, a, &lda, pivot, &info);
	}

	inline void GETRF(IntegerForLapackType ma, IntegerForLapackType na, std::complex<float>* a, IntegerForLapackType lda, IntegerForLapackType* pivot, int& info)
	{
		cgetrf_(&ma, &na, a, &lda, pivot, &info);
	}

	inline void GETRI(IntegerForLapackType na, double* a, IntegerForLapackType lda, IntegerForLapackType* pivot, double* work, IntegerForLapackType lwork, int& info)
	{
		dgetri_(&na, a, &lda, pivot, work, &lwork, &info);
	}

	inline void GETRI(IntegerForLapackType na, std::complex<double>* a, IntegerForLapackType lda, IntegerForLapackType* pivot, std::complex<double>* work, IntegerForLapackType lwork, int& info)
	{
		zgetri_(&na, a, &lda, pivot, work, &lwork, &info);
	}

	inline void GETRI(IntegerForLapackType na, float* a, IntegerForLapackType lda, IntegerForLapackType* pivot, float* work, IntegerForLapackType lwork, int& info)
	{
		sgetri_(&na, a, &lda, pivot, work, &lwork, &info);
	}

	inline void GETRI(IntegerForLapackType na, std::complex<float>* a, IntegerForLapackType lda, IntegerForLapackType* pivot, std::complex<float>* work, IntegerForLapackType lwork, int& info)
	{
		cgetri_(&na, a, &lda, pivot, work, &lwork, &info);
	}

	inline void GESDD(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, double* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  double* s,
	                  double* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  double* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  double* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  double*,
	                  // nothing
	                  IntegerForLapackType* iwork,

	                  IntegerForLapackType* info)
	{
		dgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
	}

	inline void GESDD(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, float* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  float* s,
	                  float* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  float* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  float* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  float*,
	                  // nothing
	                  IntegerForLapackType* iwork,

	                  IntegerForLapackType* info)
	{
		sgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info);
	}

	inline void GESDD(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<double>* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  double* s,
	                  std::complex<double>* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  std::complex<double>* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  std::complex<double>* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  double* rwork,
	                  // nothing
	                  IntegerForLapackType* iwork,
	                  IntegerForLapackType* info)
	{
		zgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
	}

	inline void GESDD(char* jobz, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<float>* a,
	                  // T*,
	                  IntegerForLapackType* lda,

	                  float* s,
	                  std::complex<float>* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  std::complex<float>* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  std::complex<float>* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  float* rwork,
	                  // nothing
	                  IntegerForLapackType* iwork,
	                  IntegerForLapackType* info)
	{
		cgesdd_(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, iwork, info);
	}

	inline void GESVD(char* jobz, char* jobvt, IntegerForLapackType* m, IntegerForLapackType* n, double* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  double* s,
	                  double* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  double* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  double* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  double*,
	                  // nothing
	                  IntegerForLapackType* info)
	{
		dgesvd_(jobz, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
	}

	inline void GESVD(char* jobz, char* jobvt, IntegerForLapackType* m, IntegerForLapackType* n, float* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  float* s,
	                  float* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  float* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  float* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  float*,
	                  // nothing
	                  IntegerForLapackType* info)
	{
		sgesvd_(jobz, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
	}

	inline void GESVD(char* jobz, char* jobvt, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<double>* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  double* s,
	                  std::complex<double>* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  std::complex<double>* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  std::complex<double>* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  double* rwork,
	                  // nothing
	                  IntegerForLapackType* info)
	{
		zgesvd_(jobz, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
	}

	inline void GESVD(char* jobz, char* jobvt, IntegerForLapackType* m, IntegerForLapackType* n, std::complex<float>* a,
	                  // T*,
	                  IntegerForLapackType* lda,
	                  float* s,
	                  std::complex<float>* u,
	                  // T*,
	                  IntegerForLapackType* ldu,
	                  std::complex<float>* vt,
	                  // T*,
	                  IntegerForLapackType* ldvt,
	                  std::complex<float>* work,
	                  // T*,
	                  IntegerForLapackType* lwork,
	                  float* rwork,
	                  // nothing
	                  IntegerForLapackType* info)
	{
		cgesvd_(jobz, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
	}

	inline bool isThreadSafe()
	{
		IntegerForLapackType major = 0;
		IntegerForLapackType minor = 0;
		IntegerForLapackType patch = 0;
		ilaver_(&major, &minor, &patch);
		if (major < 3)
			return false;
		return (minor >= 3);
	}
} /* namespace LAPACK */
} /* namespace psimag */

#endif // PSIMAG_LAPACK_H
