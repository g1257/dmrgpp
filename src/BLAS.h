//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for BLAS                                                     *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// * Richard Thigpen, ORNL, June 2003                                         *
// ****************************************************************************

#ifndef PSIMAG_BLAS
#define PSIMAG_BLAS
#include "AllocatorCpu.h"
#include <complex>

/** \file BLAS.h
 *  \author Thomas C. Schulthess and Richard N. Thigpen
 */

/** \brief Namespace encapsulating all PsiMag tools.
 */
namespace psimag
{

/** \brief Namespace for psimag wrappers of BLAS functions
 */
namespace BLAS
{

#ifndef PSI_BLAS_64
	typedef int IntegerForBlasType;
#else
	typedef long int IntegerForBlasType;
#endif

	//===============================================================
	// MISSING STUFF (BY G.A.)
	// ==============================================================
	extern "C" double ddot_(IntegerForBlasType*, double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	// ============================================================================
	// = Level 3 BLAS             GEMM
	// ============================================================================
	extern "C" void sgemm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dgemm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void cgemm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zgemm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);

	//*****************************************************************************
	//*                           SYMM
	//*****************************************************************************

	extern "C" void ssymm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*);

	extern "C" void dsymm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*);

	extern "C" void csymm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*);
	extern "C" void zsymm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*);

	//*****************************************************************************
	//*                           HEMM
	//*****************************************************************************

	extern "C" void chemm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zhemm_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);

	// ****************************************************************************
	// *                          SYRK
	// ****************************************************************************

	extern "C" void ssyrk_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dsyrk_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void csyrk_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zsyrk_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);

	// ****************************************************************************
	// *                          HERK
	// ****************************************************************************
	extern "C" void cherk_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zherk_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                          SYR2K
	// ****************************************************************************
	extern "C" void ssyr2k_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dsyr2k_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void csyr2k_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zsyr2k_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                          HER2k
	// ****************************************************************************
	extern "C" void cher2k_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zher2k_(char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                          TRMM
	// ****************************************************************************
	extern "C" void strmm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dtrmm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ctrmm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztrmm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                          TRSM
	// ****************************************************************************
	extern "C" void strsm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dtrsm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ctrsm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztrsm_(char*, char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *    Level 2 BLAS          GEMV
	// ****************************************************************************
	extern "C" void sgemv_(char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dgemv_(char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void cgemv_(char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zgemv_(char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                          GBMV
	// ****************************************************************************
	extern "C" void sgbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dgbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void cgbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zgbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);

	// ****************************************************************************
	// *                          HEMV
	// ****************************************************************************
	extern "C" void chemv_(char*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zhemv_(char*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                         HBMV
	// ****************************************************************************
	extern "C" void chbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zhbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ****************************************************************************
	// *                         HPMV
	// ****************************************************************************
	extern "C" void chpmv_(char*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zhpmv_(char*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         SYMV
	// ******************************************************************************
	extern "C" void ssymv_(char*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dsymv_(char*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         SBMV
	// ******************************************************************************
	extern "C" void ssbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dsbmv_(char*, IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         SPMV
	// ******************************************************************************
	extern "C" void sspmv_(char*, IntegerForBlasType*, const float*, const float*, const float*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dspmv_(char*, IntegerForBlasType*, const double*, const double*, const double*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         TRMV
	// ******************************************************************************
	extern "C" void strmv_(char*, char*, char*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dtrmv_(char*, char*, char*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ctrmv_(char*, char*, char*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztrmv_(char*, char*, char*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);

	// ******************************************************************************
	// *                         TBMV
	// ******************************************************************************
	extern "C" void stbmv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dtbmv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ctbmv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztbmv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         TPMV
	// ******************************************************************************
	extern "C" void stpmv_(char*, char*, char*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dtpmv_(char*, char*, char*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void ctpmv_(char*, char*, char*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztpmv_(char*, char*, char*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         TRSV
	// ******************************************************************************
	extern "C" void strsv_(char*, char*, char*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dtrsv_(char*, char*, char*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ctrsv_(char*, char*, char*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztrsv_(char*, char*, char*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         TBSV
	// ******************************************************************************
	extern "C" void stbsv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dtbsv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ctbsv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztbsv_(char*, char*, char*, IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         TPSV
	// ******************************************************************************
	extern "C" void stpsv_(char*, char*, char*, IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dtpsv_(char*, char*, char*, IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void ctpsv_(char*, char*, char*, IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void ztpsv_(char*, char*, char*, IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         GER
	// ******************************************************************************
	extern "C" void sger_(IntegerForBlasType*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dger_(IntegerForBlasType*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         GERU
	// ******************************************************************************
	extern "C" void cgeru_(IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zgeru_(IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         GERC
	// ******************************************************************************
	extern "C" void cgerc_(IntegerForBlasType*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zgerc_(IntegerForBlasType*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         HER
	// ******************************************************************************
	extern "C" void cher_(char*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zher_(char*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         HPR
	// ******************************************************************************
	extern "C" void chpr_(char*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*);

	extern "C" void zhpr_(char*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, std::complex<float>*);
	// ******************************************************************************
	// *                         HER2
	// ******************************************************************************
	extern "C" void cher2_(char*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zher2_(char*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         HPR2
	// ******************************************************************************
	extern "C" void chpr2_(char*, IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*);

	extern "C" void zhpr2_(char*, IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*);
	// ******************************************************************************
	// *                         SYR
	// ******************************************************************************
	extern "C" void ssyr_(char*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dsyr_(char*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         SPR
	// ******************************************************************************
	extern "C" void sspr_(char*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, float*);

	extern "C" void dspr_(char*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, double*);
	// ******************************************************************************
	// *                         SYR2
	// ******************************************************************************
	extern "C" void ssyr2_(char*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dsyr2_(char*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);
	// ******************************************************************************
	// *                         SPR2
	// ******************************************************************************
	extern "C" void sspr2_(char*, IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, const float*, IntegerForBlasType*, float*);

	extern "C" void dspr2_(char*, IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, const double*, IntegerForBlasType*, double*);
	// ******************************************************************************
	// *Level 1 BLAS
	// ******************************************************************************

	extern "C" void srotg_(float*, float*, float*, float*);
	extern "C" void drotg_(double*, double*, double*, double*);

	extern "C" void srotmg_(float*, float*, float*, float*, float*);
	extern "C" void drotmg_(double*, double*, double*, double*, double*);

	extern "C" void srot_(IntegerForBlasType*, float*, IntegerForBlasType*, float*, IntegerForBlasType*, const float*, const float*);
	extern "C" void drot_(IntegerForBlasType*, double*, IntegerForBlasType*, double*, IntegerForBlasType*, const double*, const double*);

	extern "C" void srotm_(IntegerForBlasType*, float*, IntegerForBlasType*, float*, IntegerForBlasType*, const float*);
	extern "C" void drotm_(IntegerForBlasType*, double*, IntegerForBlasType*, double*, IntegerForBlasType*, const double*);

	extern "C" void sswap_(IntegerForBlasType*, float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dswap_(IntegerForBlasType*, double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void cswap_(IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zswap_(IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);

	extern "C" void saxpy_(IntegerForBlasType*, const float*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void daxpy_(IntegerForBlasType*, const double*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void caxpy_(IntegerForBlasType*, const std::complex<float>*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zaxpy_(IntegerForBlasType*, const std::complex<double>*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);

	extern "C" void scopy_(IntegerForBlasType*, const float*, IntegerForBlasType*, float*, IntegerForBlasType*);

	extern "C" void dcopy_(IntegerForBlasType*, const double*, IntegerForBlasType*, double*, IntegerForBlasType*);

	extern "C" void ccopy_(IntegerForBlasType*, const std::complex<float>*, IntegerForBlasType*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zcopy_(IntegerForBlasType*, const std::complex<double>*, IntegerForBlasType*, std::complex<double>*, IntegerForBlasType*);

	extern "C" void sscal_(IntegerForBlasType*, const float*, float*, IntegerForBlasType*);

	extern "C" void dscal_(IntegerForBlasType*, const double*, double*, IntegerForBlasType*);

	extern "C" void cscal_(IntegerForBlasType*, const std::complex<float>*, std::complex<float>*, IntegerForBlasType*);

	extern "C" void zscal_(IntegerForBlasType*, const std::complex<double>*, std::complex<double>*, IntegerForBlasType*);

	// ============================================================================
	inline double DOT(IntegerForBlasType n, double* dx, IntegerForBlasType incx, double* dy, IntegerForBlasType incy)
	{
		return ddot_(&n, dx, &incx, dy, &incy);
	}

	// ============================================================================
	inline void GEMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, IntegerForBlasType sZ, const float& a, const float* x, IntegerForBlasType sx, const float* y, IntegerForBlasType sy, const float& b, float* z, IntegerForBlasType sz)
	{
		sgemm_(&c1, &c2, &sX, &sY, &sZ, &a, x, &sx, y, &sy, &b, z, &sz);
	}

	inline void GEMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, IntegerForBlasType sZ, const double& a, const double* x, IntegerForBlasType sx, const double* y, IntegerForBlasType sy, const double& b, double* z, IntegerForBlasType sz)
	{
		dgemm_(&c1, &c2, &sX, &sY, &sZ, &a, x, &sx, y, &sy, &b, z, &sz);
	}

	inline void GEMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, IntegerForBlasType sZ, const std::complex<float>& a, const std::complex<float>* x, IntegerForBlasType sx, const std::complex<float>* y, IntegerForBlasType sy, const std::complex<float>& b, std::complex<float>* z, IntegerForBlasType sz)
	{
		cgemm_(&c1, &c2, &sX, &sY, &sZ, &a, x, &sx, y, &sy, &b, z, &sz);
	}

	inline void GEMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, IntegerForBlasType sZ, const std::complex<double>& a, const std::complex<double>* x, IntegerForBlasType sx, const std::complex<double>* y, IntegerForBlasType sy, const std::complex<double>& b, std::complex<double>* z, IntegerForBlasType sz)
	{
		/* When  TRANSA = 'N' or 'n' then
			LDA must be at least  max( 1, m ), otherwise  LDA must be at
		least  max( 1, k ).*/

		if (c1 == 'N' || c1 == 'n') {
			if (sx < std::max(1, sX)) {
				throw PsimagLite::RuntimeError(
				    "GEMM lda < max(1, m)\n");
			}
		} else {
			if (sx < std::max(1, sZ)) {
				throw PsimagLite::RuntimeError(
				    "GEMM lda < max(1, k)\n");
			}
		}

		zgemm_(&c1, &c2, &sX, &sY, &sZ, &a, x, &sx, y, &sy, &b, z, &sz);
	}

	// ***************************************************************************
	inline void SYMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, const float& a, const float* x, IntegerForBlasType sx, const float* y, IntegerForBlasType sy, const float& b, float* z, IntegerForBlasType sz)
	{
		ssymm_(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}

	inline void SYMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, const double& a, const double* x, IntegerForBlasType sx, const double* y, IntegerForBlasType sy, const double& b, double* z, IntegerForBlasType sz)
	{
		dsymm_(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}
	inline void SYMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, const std::complex<float>& a, const std::complex<float>* x, IntegerForBlasType sx, const std::complex<float>* y, IntegerForBlasType sy, const std::complex<float>& b, std::complex<float>* z, IntegerForBlasType sz)
	{
		csymm_(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}
	inline void SYMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, const std::complex<double>& a, const std::complex<double>* x, IntegerForBlasType sx, const std::complex<double>* y, IntegerForBlasType sy, const std::complex<double>& b, std::complex<double>* z, IntegerForBlasType sz)
	{
		zsymm_(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sy, &b, z, &sz);
	}
	// ---------------------------------------------------------------------------
	inline void HEMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, const std::complex<float>& a, const std::complex<float>* x, IntegerForBlasType sx, const std::complex<float>* y, IntegerForBlasType, const std::complex<float>& b, std::complex<float>* z, IntegerForBlasType sz)
	{
		chemm_(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sx, &b, z, &sz);
	}
	inline void HEMM(char c1, char c2, IntegerForBlasType sX, IntegerForBlasType sY, const std::complex<double>& a, const std::complex<double>* x, IntegerForBlasType sx, const std::complex<double>* y, IntegerForBlasType, const std::complex<double>& b, std::complex<double>* z, IntegerForBlasType sz)
	{
		zhemm_(&c1, &c2, &sX, &sY, &a, x, &sx, y, &sx, &b, z, &sz);
	}
	// **************************************************************************
	inline void SYRK(char UPLO, char TRANS, IntegerForBlasType N, IntegerForBlasType K, const float& ALPHA, const float* A, IntegerForBlasType LDA, const float& BETA, float* C, IntegerForBlasType LDC)
	{
		ssyrk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void SYRK(char UPLO, char TRANS, IntegerForBlasType N, IntegerForBlasType K, const double& ALPHA, const double* A, IntegerForBlasType LDA, const double& BETA, double* C, IntegerForBlasType LDC)
	{
		dsyrk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void SYRK(char UPLO, char TRANS, IntegerForBlasType N, IntegerForBlasType K, const std::complex<float>& ALPHA, const std::complex<float>* A, IntegerForBlasType LDA, const std::complex<float>& BETA, std::complex<float>* C, IntegerForBlasType LDC)
	{
		csyrk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void SYRK(char UPLO, char TRANS, IntegerForBlasType N, IntegerForBlasType K, const std::complex<double>& ALPHA, const std::complex<double>* A, IntegerForBlasType LDA, const std::complex<double>& BETA, std::complex<double>* C, IntegerForBlasType LDC)
	{
		zsyrk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}

	// ***************************************************************************
	inline void HERK(char UPLO, char TRANS, IntegerForBlasType N, IntegerForBlasType K, const std::complex<float>& ALPHA, const std::complex<float>* A, IntegerForBlasType LDA, const std::complex<float>& BETA, std::complex<float>* C, IntegerForBlasType LDC)
	{
		cherk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	inline void HERK(char UPLO, char TRANS, IntegerForBlasType N, IntegerForBlasType K, const std::complex<double>& ALPHA, const std::complex<double>* A, IntegerForBlasType LDA, const std::complex<double>& BETA, std::complex<double>* C, IntegerForBlasType LDC)
	{
		zherk_(&UPLO, &TRANS, &N, &K, &ALPHA, A, &LDA, &BETA, C, &LDC);
	}
	// ***************************************************************************
	inline void SYR2K(char uplo, char trans, IntegerForBlasType n, IntegerForBlasType k, const float& alpha, const float* A, IntegerForBlasType lda, const float* B, IntegerForBlasType ldb, const float& beta, float* C, IntegerForBlasType ldc)
	{
		ssyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void SYR2K(char uplo, char trans, IntegerForBlasType n, IntegerForBlasType k, const double& alpha, const double* A, IntegerForBlasType lda, const double* B, IntegerForBlasType ldb, const double& beta, double* C, IntegerForBlasType ldc)
	{
		dsyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void SYR2k(char uplo, char trans, IntegerForBlasType n, IntegerForBlasType k, const std::complex<float>& alpha, const std::complex<float>* A, IntegerForBlasType lda, const std::complex<float>* B, IntegerForBlasType ldb, const std::complex<float>& beta, std::complex<float>* C, IntegerForBlasType ldc)
	{
		csyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void SYR2k(char uplo, char trans, IntegerForBlasType n, IntegerForBlasType k, const std::complex<double>& alpha, const std::complex<double>* A, IntegerForBlasType lda, const std::complex<double>* B, IntegerForBlasType ldb, const std::complex<double>& beta, std::complex<double>* C, IntegerForBlasType ldc)
	{
		zsyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	// ***************************************************************************
	inline void HER2k(char uplo, char trans, IntegerForBlasType n, IntegerForBlasType k, const std::complex<float>& alpha, const std::complex<float>* A, IntegerForBlasType lda, const std::complex<float>* B, IntegerForBlasType ldb, const std::complex<float>& beta, std::complex<float>* C, IntegerForBlasType ldc)
	{
		cher2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	inline void HER2k(char uplo, char trans, IntegerForBlasType n, IntegerForBlasType k, const std::complex<double>& alpha, const std::complex<double>* A, IntegerForBlasType lda, const std::complex<double>* B, IntegerForBlasType ldb, const std::complex<double>& beta, std::complex<double>* C, IntegerForBlasType ldc)
	{
		zher2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
	}
	// ********************************************************************************
	inline void TRMM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const float& alpha, const float* A, IntegerForBlasType lda, float* B, IntegerForBlasType ldb)
	{
		strmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void TRMM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const double& alpha, const double* A, IntegerForBlasType lda, double* B, IntegerForBlasType ldb)
	{
		dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void TRMM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* A, IntegerForBlasType lda, std::complex<float>* B, IntegerForBlasType ldb)
	{
		ctrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void TRMM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* A, IntegerForBlasType lda, std::complex<double>* B, IntegerForBlasType ldb)
	{
		ztrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	// ********************************************************************************
	inline void TRSM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const float& alpha, const float* A, IntegerForBlasType lda, float* B, IntegerForBlasType ldb)
	{
		strsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void TRSM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const double& alpha, const double* A, IntegerForBlasType lda, double* B, IntegerForBlasType ldb)
	{
		dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void TRSM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* A, IntegerForBlasType lda, std::complex<float>* B, IntegerForBlasType ldb)
	{
		ctrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	inline void TRSM(char side, char uplo, char transa, char diag, IntegerForBlasType m, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* A, IntegerForBlasType lda, std::complex<double>* B, IntegerForBlasType ldb)
	{
		ztrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb);
	}
	// ***************************************************************************

	inline void GEMV(char c, IntegerForBlasType M, IntegerForBlasType N, const float& alpha, const float* A, IntegerForBlasType ldA, const float* x, IntegerForBlasType incX, const float& beta, float* y, IntegerForBlasType incY)
	{
		sgemv_(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ----------------------------------------------------------------------------
	inline void GEMV(char c, IntegerForBlasType M, IntegerForBlasType N, const double& alpha, const double* A, IntegerForBlasType ldA, const double* x, IntegerForBlasType incX, const double& beta, double* y, IntegerForBlasType incY)
	{
		dgemv_(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ---------------------------------------------------------------------------
	inline void GEMV(char c, IntegerForBlasType M, IntegerForBlasType N, const std::complex<float>& alpha, const std::complex<float>* A, IntegerForBlasType ldA, const std::complex<float>* x, IntegerForBlasType incX, const std::complex<float>& beta, std::complex<float>* y, IntegerForBlasType incY)
	{
		cgemv_(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ---------------------------------------------------------------------------
	inline void GEMV(char c, IntegerForBlasType M, IntegerForBlasType N, const std::complex<double>& alpha, const std::complex<double>* A, IntegerForBlasType ldA, const std::complex<double>* x, IntegerForBlasType incX, const std::complex<double>& beta, std::complex<double>* y, IntegerForBlasType incY)
	{
		zgemv_(&c, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
	}
	// ----------------------------------------------------------------------------
	inline void GBMV(char trans, IntegerForBlasType m, IntegerForBlasType n, IntegerForBlasType kl, IntegerForBlasType ku, const float& alpha, const float* A, IntegerForBlasType lda, const float* x, IntegerForBlasType incx, const float& beta, float* y, IntegerForBlasType incy)
	{
		sgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	inline void GBMV(char trans, IntegerForBlasType m, IntegerForBlasType n, IntegerForBlasType kl, IntegerForBlasType ku, const double& alpha, const double* A, IntegerForBlasType lda, const double* x, IntegerForBlasType incx, const double& beta, double* y, IntegerForBlasType incy)
	{
		dgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	inline void GBMV(char trans, IntegerForBlasType m, IntegerForBlasType n, IntegerForBlasType kl, IntegerForBlasType ku, const std::complex<float>& alpha, const std::complex<float>* A, IntegerForBlasType lda, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>& beta, std::complex<float>* y, IntegerForBlasType incy)
	{
		cgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	inline void GBMV(char trans, IntegerForBlasType m, IntegerForBlasType n, IntegerForBlasType kl, IntegerForBlasType ku, const std::complex<double>& alpha, const std::complex<double>* A, IntegerForBlasType lda, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>& beta, std::complex<double>* y, IntegerForBlasType incy)
	{
		zgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void HEMV(char uplo, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* a, IntegerForBlasType lda, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>& beta, std::complex<float>* y, IntegerForBlasType incy)
	{
		chemv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void HEMV(char uplo, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* a, IntegerForBlasType lda, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>& beta, std::complex<double>* y, IntegerForBlasType incy)
	{
		zhemv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// **************************************************************************
	inline void HBMV(char uplo, IntegerForBlasType n, IntegerForBlasType k, const std::complex<float>& alpha, const std::complex<float>* a, IntegerForBlasType lda, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>& beta, std::complex<float>* y, IntegerForBlasType incy)
	{
		chbmv_(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void HBMV(char uplo, IntegerForBlasType n, IntegerForBlasType k, const std::complex<double>& alpha, const std::complex<double>* a, IntegerForBlasType lda, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>& beta, std::complex<double>* y, IntegerForBlasType incy)
	{
		zhbmv_(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// ***************************************************************************
	inline void HPMV(char uplo, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* ap, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>& beta, std::complex<float>* y, IntegerForBlasType incy)
	{
		chpmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	inline void HPMV(char uplo, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* ap, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>& beta, std::complex<double>* y, IntegerForBlasType incy)
	{
		zhpmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	// ***************************************************************************
	inline void SYMV(char uplo, IntegerForBlasType n, const float& alpha, const float* a, IntegerForBlasType lda, const float* x, IntegerForBlasType incx, const float& beta, float* y, IntegerForBlasType incy)
	{
		ssymv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void SYMV(char uplo, IntegerForBlasType n, const double& alpha, const double* a, IntegerForBlasType lda, const double* x, IntegerForBlasType incx, const double& beta, double* y, IntegerForBlasType incy)
	{
		dsymv_(&uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void SBMV(char uplo, IntegerForBlasType n, IntegerForBlasType k, const float& alpha, const float* a, IntegerForBlasType lda, const float* x, IntegerForBlasType incx, const float& beta, float* y, IntegerForBlasType incy)
	{
		ssbmv_(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	inline void SBMV(char uplo, IntegerForBlasType n, IntegerForBlasType k, const double& alpha, const double* a, IntegerForBlasType lda, const double* x, IntegerForBlasType incx, const double& beta, double* y, IntegerForBlasType incy)
	{
		dsbmv_(&uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void SPMV(char uplo, IntegerForBlasType n, const float& alpha, const float* ap, const float* x, IntegerForBlasType incx, const float& beta, float* y, IntegerForBlasType incy)
	{
		sspmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	inline void SPMV(char uplo, IntegerForBlasType n, const double& alpha, const double* ap, const double* x, IntegerForBlasType incx, const double& beta, double* y, IntegerForBlasType incy)
	{
		dspmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy);
	}
	// ****************************************************************************
	inline void TRMV(char uplo, char trans, char diag, IntegerForBlasType n, const float* a, IntegerForBlasType lda, float* x, IntegerForBlasType incx)
	{
		strmv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void TRMV(char uplo, char trans, char diag, IntegerForBlasType n, const double* a, IntegerForBlasType lda, double* x, IntegerForBlasType incx)
	{
		dtrmv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void TRMV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<float>* a, IntegerForBlasType lda, std::complex<float>* x, IntegerForBlasType incx)
	{
		ctrmv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void TRMV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<double>* a, IntegerForBlasType lda, std::complex<double>* x, IntegerForBlasType incx)
	{
		ztrmv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void TBMV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const float* a, IntegerForBlasType lda, float* x, IntegerForBlasType incx)
	{
		stbmv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void TBMV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const double* a, IntegerForBlasType lda, double* x, IntegerForBlasType incx)
	{
		dtbmv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void TBMV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const std::complex<float>* a, IntegerForBlasType lda, std::complex<float>* x, IntegerForBlasType incx)
	{
		ctbmv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void TBMV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const std::complex<double>* a, IntegerForBlasType lda, std::complex<double>* x, IntegerForBlasType incx)
	{
		ztbmv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void TPMV(char uplo, char trans, char diag, IntegerForBlasType n, const float* ap, float* x, IntegerForBlasType incx)
	{
		stpmv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void TPMV(char uplo, char trans, char diag, IntegerForBlasType n, const double* ap, double* x, IntegerForBlasType incx)
	{
		dtpmv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void TPMV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<float>* ap, std::complex<float>* x, IntegerForBlasType incx)
	{
		ctpmv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void TPMV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<double>* ap, std::complex<double>* x, IntegerForBlasType incx)
	{
		ztpmv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	// ****************************************************************************
	inline void TRSV(char uplo, char trans, char diag, IntegerForBlasType n, const float* a, IntegerForBlasType lda, float* x, IntegerForBlasType incx)
	{
		strsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void TRSV(char uplo, char trans, char diag, IntegerForBlasType n, const double* a, IntegerForBlasType lda, double* x, IntegerForBlasType incx)
	{
		dtrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void TRSV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<float>* a, IntegerForBlasType lda, std::complex<float>* x, IntegerForBlasType incx)
	{
		ctrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	inline void TRSV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<double>* a, IntegerForBlasType lda, std::complex<double>* x, IntegerForBlasType incx)
	{
		ztrsv_(&uplo, &trans, &diag, &n, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void TBSV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const float* a, IntegerForBlasType lda, float* x, IntegerForBlasType incx)
	{
		stbsv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void TBSV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const double* a, IntegerForBlasType lda, double* x, IntegerForBlasType incx)
	{
		dtbsv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void TBSV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const std::complex<float>* a, IntegerForBlasType lda, std::complex<float>* x, IntegerForBlasType incx)
	{
		ctbsv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	inline void TBSV(char uplo, char trans, char diag, IntegerForBlasType n, IntegerForBlasType k, const std::complex<double>* a, IntegerForBlasType lda, std::complex<double>* x, IntegerForBlasType incx)
	{
		ztbsv_(&uplo, &trans, &diag, &n, &k, a, &lda, x, &incx);
	}
	// ****************************************************************************
	inline void TPSV(char uplo, char trans, char diag, IntegerForBlasType n, const float* ap, float* x, IntegerForBlasType incx)
	{
		stpsv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void TPSV(char uplo, char trans, char diag, IntegerForBlasType n, const double* ap, double* x, IntegerForBlasType incx)
	{
		dtpsv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void TPSV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<float>* ap, std::complex<float>* x, IntegerForBlasType incx)
	{
		ctpsv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	inline void TPSV(char uplo, char trans, char diag, IntegerForBlasType n, const std::complex<double>* ap, std::complex<double>* x, IntegerForBlasType incx)
	{
		ztpsv_(&uplo, &trans, &diag, &n, ap, x, &incx);
	}
	// ****************************************************************************
	inline void GER(IntegerForBlasType m, IntegerForBlasType n, const float& alpha, const float* x, IntegerForBlasType incx, const float* y, IntegerForBlasType incy, float* a, IntegerForBlasType lda)
	{
		sger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void GER(IntegerForBlasType m, IntegerForBlasType n, const double& alpha, const double* x, IntegerForBlasType incx, const double* y, IntegerForBlasType incy, double* a, IntegerForBlasType lda)
	{
		dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// ****************************************************************************
	inline void GERU(IntegerForBlasType m, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>* y, IntegerForBlasType incy, std::complex<float>* a, IntegerForBlasType lda)
	{
		cgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void GERU(IntegerForBlasType m, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>* y, IntegerForBlasType incy, std::complex<double>* a, IntegerForBlasType lda)
	{
		zgeru_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// ****************************************************************************
	inline void GERC(IntegerForBlasType m, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>* y, IntegerForBlasType incy, std::complex<float>* a, IntegerForBlasType lda)
	{
		cgerc_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void GERC(IntegerForBlasType m, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>* y, IntegerForBlasType incy, std::complex<double>* a, IntegerForBlasType lda)
	{
		zgerc_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// *****************************************************************************
	inline void HER(char uplo, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* x, IntegerForBlasType incx, std::complex<float>* a, IntegerForBlasType lda)
	{
		cher_(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	inline void HER(char uplo, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* x, IntegerForBlasType incx, std::complex<double>* a, IntegerForBlasType lda)
	{
		zher_(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	// *****************************************************************************
	inline void HPR(char uplo, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* x, IntegerForBlasType incx, std::complex<float>* ap)
	{
		chpr_(&uplo, &n, &alpha, x, &incx, ap);
	}
	inline void HPR(char uplo, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* x, IntegerForBlasType incx, std::complex<float>* ap)
	{
		zhpr_(&uplo, &n, &alpha, x, &incx, ap);
	}
	// *****************************************************************************
	inline void HER2(char uplo, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>* y, IntegerForBlasType incy, std::complex<float>* a, IntegerForBlasType lda)
	{
		cher2_(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void HER2(char uplo, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>* y, IntegerForBlasType incy, std::complex<double>* a, IntegerForBlasType lda)
	{
		zher2_(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// *****************************************************************************
	inline void HPR2(char uplo, IntegerForBlasType n, const std::complex<float>& alpha, const std::complex<float>* x, IntegerForBlasType incx, const std::complex<float>* y, IntegerForBlasType incy, std::complex<float>* ap)
	{
		chpr2_(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}
	inline void HPR2(char uplo, IntegerForBlasType n, const std::complex<double>& alpha, const std::complex<double>* x, IntegerForBlasType incx, const std::complex<double>* y, IntegerForBlasType incy, std::complex<double>* ap)
	{
		zhpr2_(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}
	// *****************************************************************************
	inline void SYR(char uplo, IntegerForBlasType n, const float& alpha, const float* x, IntegerForBlasType incx, float* a, IntegerForBlasType lda)
	{
		ssyr_(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	inline void SYR(char uplo, IntegerForBlasType n, const double& alpha, const double* x, IntegerForBlasType incx, double* a, IntegerForBlasType lda)
	{
		dsyr_(&uplo, &n, &alpha, x, &incx, a, &lda);
	}
	// ****************************************************************************
	inline void SPR(char uplo, IntegerForBlasType n, const float& alpha, const float* x, IntegerForBlasType incx, float* ap)
	{
		sspr_(&uplo, &n, &alpha, x, &incx, ap);
	}
	inline void SPR(char uplo, IntegerForBlasType n, const double& alpha, const double* x, IntegerForBlasType incx, double* ap)
	{
		dspr_(&uplo, &n, &alpha, x, &incx, ap);
	}
	// ****************************************************************************
	inline void SYR2(char uplo, IntegerForBlasType n, const float& alpha, const float* x,

	    IntegerForBlasType incx,
	    const float* y,
	    IntegerForBlasType incy,
	    float* a,
	    IntegerForBlasType lda)
	{
		ssyr2_(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	inline void SYR2(char uplo, IntegerForBlasType n, const double& alpha, const double* x,

	    IntegerForBlasType incx,
	    const double* y,
	    IntegerForBlasType incy,
	    double* a,
	    IntegerForBlasType lda)
	{
		dsyr2_(&uplo, &n, &alpha, x, &incx, y, &incy, a, &lda);
	}
	// ****************************************************************************
	inline void SPR2(char uplo, IntegerForBlasType n, const float& alpha, const float* x,

	    IntegerForBlasType incx,
	    const float* y,
	    IntegerForBlasType incy,
	    float* ap)
	{
		sspr2_(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}
	inline void SPR2(char uplo, IntegerForBlasType n, const double& alpha, const double* x,

	    IntegerForBlasType incx,
	    const double* y,
	    IntegerForBlasType incy,
	    double* ap)
	{
		dspr2_(&uplo, &n, &alpha, x, &incx, y, &incy, ap);
	}

	// ****************************************************************************

	inline void AXPY(IntegerForBlasType size, const float& a, const float* x, IntegerForBlasType sx, float* y, IntegerForBlasType sy)
	{
		saxpy_(&size, &a, x, &sx, y, &sy);
	}

	inline void AXPY(IntegerForBlasType size, const double& a, const double* x, IntegerForBlasType sx, double* y, IntegerForBlasType sy)
	{
		daxpy_(&size, &a, x, &sx, y, &sy);
	}

	inline void AXPY(IntegerForBlasType size, const std::complex<float>& a, const std::complex<float>* x, IntegerForBlasType sx,

	    std::complex<float>* y,
	    IntegerForBlasType sy)
	{
		caxpy_(&size, &a, x, &sx, y, &sy);
	}

	inline void AXPY(IntegerForBlasType size, const std::complex<double>& a, const std::complex<double>* x, IntegerForBlasType sx,

	    std::complex<double>* y,
	    IntegerForBlasType sy)
	{
		zaxpy_(&size, &a, x, &sx, y, &sy);
	}
	// ----------------------------------------------------------------------------
	inline void COPY(IntegerForBlasType size, const float* x, IntegerForBlasType sx, float* y, IntegerForBlasType sy)
	{
		scopy_(&size, x, &sx, y, &sy);
	}

	inline void COPY(IntegerForBlasType size, const double* x, IntegerForBlasType sx, double* y, IntegerForBlasType sy)
	{
		dcopy_(&size, x, &sx, y, &sy);
	}

	inline void COPY(IntegerForBlasType size, const std::complex<float>* x, IntegerForBlasType sx, std::complex<float>* y, IntegerForBlasType sy)
	{
		ccopy_(&size, x, &sx, y, &sy);
	}

	inline void COPY(IntegerForBlasType size, const std::complex<double>* x, IntegerForBlasType sx, std::complex<double>* y, IntegerForBlasType sy)
	{
		zcopy_(&size, x, &sx, y, &sy);
	}
	// ----------------------------------------------------------------------------
	inline void SCAL(IntegerForBlasType size, const float& a, float* y, IntegerForBlasType sy)
	{
		sscal_(&size, &a, y, &sy);
	}

	inline void SCAL(IntegerForBlasType size, const double& a, double* y, IntegerForBlasType sy)
	{
		dscal_(&size, &a, y, &sy);
	}

	inline void SCAL(IntegerForBlasType size, const std::complex<float>& a, std::complex<float>* y, IntegerForBlasType sy)
	{
		cscal_(&size, &a, y, &sy);
	}

	inline void SCAL(IntegerForBlasType size, const std::complex<double>& a, std::complex<double>* y, IntegerForBlasType sy)
	{
		zscal_(&size, &a, y, &sy);
	}
} /* namespace BLAS */
} /* namespace psimag */

#endif /* PSIMAG_BLAS */
