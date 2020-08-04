#ifndef INTER_NODE_H
#define INTER_NODE_H
#include "Vector.h"
#include "Concurrency.h"

#ifdef USE_MPI
#include "InterNodeMpi.h"
#else
#include "InterNodeSerial.h"
#endif

#endif // INTER_NODE_H
