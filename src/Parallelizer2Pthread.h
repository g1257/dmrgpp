#ifndef PARALLELIZER2PTHREAD_H
#define PARALLELIZER2PTHREAD_H
#include <pthread.h>
#include <iostream>
#include <algorithm>
#include "Vector.h"
#include <sched.h>
#include <unistd.h>
#include "TypeToString.h"
#include "CodeSectionParams.h"
#include "LoadBalancerDefault.h"

#ifdef __linux__
#include <sys/types.h>
#include <sys/syscall.h>
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifdef _GNU_SOURCE
#include <errno.h>
#include <string.h>
#endif

template<typename SomeLambdaType,
         typename LoadBalancerType=PsimagLite::LoadBalancerDefault
         >
struct PthreadFunctionStruct {
	PthreadFunctionStruct()
	    : pfh(0),loadBalancer(0),threadNum(0),nthreads(0),start(0),end(0),cpu(0)
	{}

	const SomeLambdaType* pfh;
	const LoadBalancerType* loadBalancer;
	int threadNum;
	SizeType nthreads;
	SizeType start;
	SizeType end;
	SizeType cpu;
};

template<typename SomeLambdaType>
void *thread_function_wrapper(void *dummyPtr)
{
	PthreadFunctionStruct<SomeLambdaType> *pfs =
	        static_cast<PthreadFunctionStruct<SomeLambdaType> *>(dummyPtr);

	const SomeLambdaType* pfh = pfs->pfh;

	int s = 0;
#ifdef __linux__
	s = sched_getcpu();
#endif
	if (s >= 0) pfs->cpu = s;

	SizeType blockSize = pfs->loadBalancer->blockSize(pfs->threadNum);

	for (SizeType p = 0; p < blockSize; ++p) {
		SizeType taskNumber = pfs->loadBalancer->taskNumber(pfs->threadNum, p);
		if (taskNumber + pfs->start >= pfs->end) break;
		(*pfh)(taskNumber + pfs->start, pfs->threadNum);
	}

	return 0;
}

namespace PsimagLite {

template<typename LoadBalancerType=LoadBalancerDefault>
class Parallizer2 {

public:

	typedef LoadBalancerDefault::VectorSizeType VectorSizeType;

	Parallizer2(SizeType nthreads)
	    : nthreads_(nthreads) {}

	// no weights, no balancer ==> create weights, set all weigths to 1, delegate
	template<typename SomeLambdaType>
	void parallelFor(SizeType start, SizeType end, const SomeLambdaType& lambda)
	{
		VectorSizeType weights(end - start, 1);
		parallelFor(start, end, lambda, weights);
	}

	// weights, no balancer ==> create balancer with weights ==> delegate
	template<typename SomeLambdaType>
	void parallelFor(SizeType start,
	                 SizeType end,
	                 const SomeLambdaType& lambda,
	                 const VectorSizeType& weights)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights, nthreads_);
		parallelFor(start, end, lambda, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	template<typename SomeLambdaType>
	void parallelFor(SizeType start,
	                 SizeType end,
	                 const SomeLambdaType& lambda,
	                 const LoadBalancerType& loadBalancer)
	{
		PthreadFunctionStruct<SomeLambdaType>* pfs =
		new PthreadFunctionStruct<SomeLambdaType>[nthreads_];
		pthread_t* thread_id = new pthread_t[nthreads_];
		pthread_attr_t** attr = new pthread_attr_t*[nthreads_];

		for (SizeType j = 0; j < nthreads_; ++j) {
			pfs[j].pfh = &lambda;
			pfs[j].loadBalancer = &loadBalancer;
			pfs[j].threadNum = j;
			pfs[j].start = start;
			pfs[j].end = end;
			pfs[j].nthreads = nthreads_;

			attr[j] = new pthread_attr_t;

			int ret = pthread_attr_init(attr[j]);
			checkForError(ret);

			ret = pthread_create(&thread_id[j],
			                     attr[j],
			                     thread_function_wrapper<SomeLambdaType>,
			                     &pfs[j]);
			checkForError(ret);
		}

		for (SizeType j=0; j <nthreads_; ++j) pthread_join(thread_id[j], 0);
		for (SizeType j=0; j <nthreads_; ++j) {
			int ret = pthread_attr_destroy(attr[j]);
			checkForError(ret);
			delete attr[j];
			attr[j] = 0;
		}

		delete [] attr;
		delete [] thread_id;
		delete [] pfs;
	}

private:

	void checkForError(int ret) const
	{
		if (ret == 0) return;
		std::cerr<<"PthreadsNg ERROR: "<<strerror(ret)<<"\n";
	}

	SizeType nthreads_;
};
}
#endif // PARALLELIZER2PTHREAD_H
