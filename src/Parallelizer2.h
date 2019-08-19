#ifndef PARALLELIZER2_H
#define PARALLELIZER2_H
#include "Vector.h"

#ifdef USE_PTHREADS
#include "Parallelizer2Pthread.h"
#else
#include "Parallelizer2Serial.h"
#endif

#endif // PARALLELIZER2_H
