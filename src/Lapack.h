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

extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
        std::complex<double> *,int *, double *, int *);
extern "C" void dsyev_(char *,char *,int *,double *,int *, double *,double *,int *,int *);

#endif
