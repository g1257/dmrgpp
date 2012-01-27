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
  extern "C" int  dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);

  extern "C" int  zgetrf_(int* M, int* N, std::complex<double>* A, int* LDA, int* IPIV, int* INFO);

  extern "C" int  zgetri_(int* N, std::complex<double>* A, int* LDA, int* IPIV,  std::complex<double>* WORK, int* LWORK, int* INFO);

  extern "C" int  dgetri_(int* N, double* A, int* LDA, int* IPIV,  double* WORK, int* LWORK, int* INFO);

// ============================================================================
  inline void GESV(int ma,int mb,float* a,int lda,int* pivot,
		   float* b,int ldb,int& info) {
    sgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }
  inline void GESV(int ma,int mb,double* a,int lda,int* pivot,
		   double* b,int ldb,int& info) {
    dgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }
  inline void GESV(int ma,int mb,std::complex<float>* a,int lda,int* pivot,
		   std::complex<float>* b,int ldb,int& info) {
    cgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }
  inline void GESV(int ma,int mb,std::complex<double>* a,int lda,int* pivot,
		   std::complex<double>* b,int ldb,int& info) {
    zgesv_(&ma,&mb,a,&lda,pivot,b,&ldb,&info);
  }

  inline void GETRF(int ma, int na, double* a, int lda, int* pivot, int& info) {
    dgetrf_(&ma,&na,a,&lda,pivot,&info);
  }

  inline void GETRF(int ma,int na,std::complex<double>* a,int lda,int* pivot,int& info) {
	  zgetrf_(&ma,&na,a,&lda,pivot,&info);
  }
  
  inline void GETRI(int na, double* a, int lda, int* pivot, double* work, int lwork, int& info) {
    dgetri_(&na,a,&lda,pivot,work,&lwork,&info);
  }

  inline void GETRI(int na, std::complex<double>* a, int lda, int* pivot, std::complex<double>* work, int lwork, int& info) {
    zgetri_(&na,a,&lda,pivot,work,&lwork,&info);
  }
  
}      /* namespace LAPACK */
}      /* namespace psimag */

#endif // PSIMAG_LAPACK_H
