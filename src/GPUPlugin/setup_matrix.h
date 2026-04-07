#ifndef SETUP_MATRIX_H
#define SETUP_MATRIX_H 

#include "dmrg_types.h"
#include "Vector.h"

typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

#if defined(USE_COMPLEX_Z)
std::complex<double> makeFloat(double zr,double zi);
#else
double makeFloat(double zr,double zi);
#endif

template<typename T>
void setup_matrix(SizeType noperator,
                  SizeType npatches,
		  std::vector<T>& Abatch,
                  const VectorSizeType& left_patch_size_,
                  std::vector<T*>& Amatrix_,
                  VectorSizeType& ld_Amatrix_);

#endif

