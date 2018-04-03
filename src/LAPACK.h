//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for BLAS                                                     *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// ****************************************************************************

#ifndef PSIMAG_LAPACK
#define PSIMAG_LAPACK

#include <complex>
#include "Matrix.h"

/** \file LAPACK
 *  \author Thomas C. Schulthess, MSS
 */

namespace psimag {

/** \brief Namespace encapsulating psimag wrappers of LAPACK functions.
 */
namespace LAPACK {

// ============================================================================
extern "C" void sgesv_(int*,int*,float*,int*,int*,float*,int*,int*);

extern "C" void dgesv_(int*,int*,double*,int*,int*,double*,int*,int*);

extern "C" void cgesv_(int*,int*,std::complex<float>*,int*,int*,std::complex<float>*,
                       int*,int*);
extern "C" void zgesv_(int*,int*,std::complex<double>*,int*,int*,std::complex<double>*,
                       int*,int*);

//MSS
extern "C" int  dgetrf_(int*, int*, double*, int*, int*, int*);

extern "C" int  zgetrf_(int*, int*, std::complex<double>*, int*, int*, int*);

extern "C" int  zgetri_(int*,
                        std::complex<double>*,
                        int*,
                        int*,
                        std::complex<double>*,
                        int*,
                        int*);

extern "C" int  dgetri_(int*, double*, int*, int*,  double*, int*, int*);

extern "C" int  dgesdd_(char* jobz,
                        int* m,
                        int* n,
                        double* a,
                        int* lda,
                        double* s,
                        double* u,
                        int* ldu,
                        double* vt,
                        int* ldvt,
                        double* work,
                        int* lwork,
                        int* iwork,
                        int* info);

extern "C" int  sgesdd_(char* jobz,
                        int* m,
                        int* n,
                        float* a,
                        int* lda,
                        float* s,
                        float* u,
                        int* ldu,
                        float* vt,
                        int* ldvt,
                        float* work,
                        int* lwork,
                        int* iwork,
                        int* info);

extern "C" int  zgesdd_(char* jobz,
                        int* m,
                        int* n,
                        std::complex<double>* a,
                        int* lda,
                        double* s,
                        std::complex<double>* u,
                        int* ldu,
                        std::complex<double>* vt,
                        int* ldvt,
                        std::complex<double>* work,
                        int* lwork,
                        double* rwork,
                        int* iwork,
                        int* info);

extern "C" int  cgesdd_(char* jobz,
                        int* m,
                        int* n,
                        std::complex<float>* a,
                        int* lda,
                        float* s,
                        std::complex<float>* u,
                        int* ldu,
                        std::complex<float>* vt,
                        int* ldvt,
                        std::complex<float>* work,
                        int* lwork,
                        float* rwork,
                        int* iwork,
                        int* info);

extern "C" void ilaver_(int*, int*, int*);

// ============================================================================
inline void GESV(int ma,
                 int mb,
                 float* a,
                 int lda,
                 int* pivot,
                 float* b,
                 int ldb,
                 int& info)
{
	sgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
}

inline void GESV(int ma,
                 int mb,
                 double* a,
                 int lda,
                 int* pivot,
                 double* b,
                 int ldb,
                 int& info)
{
	dgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
}

inline void GESV(int ma,
                 int mb,
                 std::complex<float>* a,
                 int lda,
                 int* pivot,
                 std::complex<float>* b,
                 int ldb,
                 int& info)
{
	cgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
}

inline void GESV(int ma,
                 int mb,
                 std::complex<double>* a,
                 int lda,
                 int* pivot,
                 std::complex<double>* b,
                 int ldb,
                 int& info)
{
	zgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
}

inline void GETRF(int ma, int na, double* a, int lda, int* pivot, int& info)
{
	dgetrf_(&ma,&na,a,&lda,pivot,&info);
}

inline void GETRF(int ma,int na,std::complex<double>* a,int lda,int* pivot,int& info)
{
	zgetrf_(&ma,&na,a,&lda,pivot,&info);
}

inline void GETRI(int na,
                  double* a,
                  int lda,
                  int* pivot,
                  double* work,
                  int lwork,
                  int& info)
{
	dgetri_(&na,a,&lda,pivot,work,&lwork,&info);
}

inline void GETRI(int na,
                  std::complex<double>* a,
                  int lda,
                  int* pivot,
                  std::complex<double>* work,
                  int lwork,
                  int& info)
{
	zgetri_(&na,a,&lda,pivot,work,&lwork,&info);
}

inline bool isThreadSafe()
{
	int major = 0;
	int minor = 0;
	int patch = 0;
	ilaver_(&major, &minor, &patch);
	if (major < 3) return false;
	return (minor >= 3);
}
}      /* namespace LAPACK */
}      /* namespace psimag */

#endif // PSIMAG_LAPACK_H
