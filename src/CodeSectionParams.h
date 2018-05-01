#ifndef CODESECTION_PARAMS_H
#define CODESECTION_PARAMS_H
#include "AllocatorCpu.h"

namespace PsimagLite {

struct CodeSectionParams {

	CodeSectionParams(SizeType threads,
	            bool a = false,
	            size_t s = 0)
	    : npthreads(threads), setAffinities(a), stackSize(s)
	{}

	SizeType npthreads;
	bool setAffinities;
	size_t stackSize;
};
}
#endif // CODESECTION_PARAMS_H

