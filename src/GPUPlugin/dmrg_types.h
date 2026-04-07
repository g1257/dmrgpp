#ifndef DMRG_TYPES_H
#define DMRG_TYPES_H

typedef int IntegerType;

#if defined(USE_COMPLEX_Z)

#include <complex.h>

typedef std::complex<double> MYTYPE;

#else

typedef double MYTYPE;

#endif

#endif
