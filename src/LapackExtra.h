// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
// END LICENSE BLOCK

#ifndef LAPACK_H_
#define LAPACK_H_
#include <complex>

#ifndef PSI_LAPACK_64
typedef int IntegerForLapackType;
#else
typedef long int IntegerForLapackType;
#endif

extern "C" void   zheev_(char*,
                         char*,
                         IntegerForLapackType*,
                         std::complex<double>*,
                         IntegerForLapackType*,
                         double*,
                         std::complex<double>*,
                         IntegerForLapackType*,
                         double*,
                         IntegerForLapackType*);
extern "C" void   cheev_(char*,
                         char*,
                         IntegerForLapackType*,
                         std::complex<float>*,
                         IntegerForLapackType*,
                         float*,
                         std::complex<float>*,
                         IntegerForLapackType*,
                         float*,
                         IntegerForLapackType*);
extern "C" void dsyev_(char*,
                       char*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*,
                       double*,
                       double*,
                       IntegerForLapackType*,
                       IntegerForLapackType*);
extern "C" void ssyev_(char*,
                       char*,
                       IntegerForLapackType*,
                       float*,
                       IntegerForLapackType*,
                       float*,
                       float*,
                       IntegerForLapackType*,
                       IntegerForLapackType*);

extern "C" void zgeev_(char*,
                       char*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       std::complex<double>*,
                       IntegerForLapackType*,
                       double*,
                       IntegerForLapackType*);
#endif
