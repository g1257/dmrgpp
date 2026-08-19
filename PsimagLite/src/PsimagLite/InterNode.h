#ifndef INTER_NODE_H
#define INTER_NODE_H
#include "Concurrency.h"
#include "Vector.h"

#ifdef DMRGPP_USE_MPI
#include "InterNodeMpi.h"
#else
#include "InterNodeSerial.h"
#endif

#endif // INTER_NODE_H
