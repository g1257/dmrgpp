#ifndef DMRG_TYPES_H
#define DMRG_TYPES_H

using IntegerType = int;

#if defined(USE_COMPLEX_Z)

#include <complex>

typedef std::complex<double> MYTYPE;

#else

typedef double MYTYPE;

#endif

#endif
