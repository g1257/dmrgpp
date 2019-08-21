#ifndef PARALLELIZER2_H
#define PARALLELIZER2_H
#include "Vector.h"

#ifdef USE_PTHREADS
#include "Parallelizer2Pthread.h"
#else

#ifdef _OPENMP
#include "Parallelizer2OpenMP.h"
#else
#include "Parallelizer2Serial.h"
#endif

#endif

#endif // PARALLELIZER2_H
