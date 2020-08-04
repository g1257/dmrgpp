#ifndef INTER_NODE_MPI_H
#define INTER_NODE_MPI_H
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include "Vector.h"
#include <sched.h>
#include <unistd.h>
#include "TypeToString.h"
#include "LoadBalancerDefault.h"

template<typename SomeLambdaType,
         typename LoadBalancerType=PsimagLite::LoadBalancerDefault
         >
struct MpiFunctionStruct {
	MpiFunctionStruct()
	    : pfh(0),loadBalancer(0),nsize(0),start(0),end(0)
	{}

	const SomeLambdaType* pfh;
	const LoadBalancerType* loadBalancer;
	SizeType nsize;
	SizeType start;
	SizeType end;
};

template<typename SomeLambdaType, typename SomeLoadBalancer>
void *mpiFunctionWrapper(void *dummyPtr)
{
	MpiFunctionStruct<SomeLambdaType, SomeLoadBalancer> *pfs =
	        static_cast<MpiFunctionStruct<SomeLambdaType, SomeLoadBalancer> *>(dummyPtr);

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

	int retval = 0;
	pthread_exit(static_cast<void*>(&retval));
	return 0;
}

namespace PsimagLite {

template<typename LoadBalancerType=LoadBalancerDefault>
class Parallelizer2 {

public:

	typedef LoadBalancerDefault::VectorSizeType VectorSizeType;

	Parallelizer2(const CodeSectionParams& codeParams)
	    : nthreads_(codeParams.npthreads), stackSize_(codeParams.stackSize) {}

	SizeType numberOfThreads() const { return nthreads_; }

	String name() const { return "pthread"; }

	// no weights, no balancer ==> create weights, set all weigths to 1, delegate
	template<typename SomeLambdaType>
	void parallelFor(SizeType start, SizeType end, const SomeLambdaType& lambda)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(end - start, nthreads_);
		parallelFor(start, end, lambda, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// weights, no balancer ==> create balancer with weights ==> delegate
	template<typename SomeLambdaType>
	void parallelFor(SizeType start,
	                 SizeType end,
	                 const SomeLambdaType& lambda,
	                 const VectorSizeType& weights)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights.size(), nthreads_);
		loadBalancer->setWeights(weights);
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
		MpiFunctionStruct<SomeLambdaType, LoadBalancerType>* pfs =
		new MpiFunctionStruct<SomeLambdaType, LoadBalancerType>[nthreads_];
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
			int ret = (stackSize_ > 0) ? pthread_attr_setstacksize(attr[j], stackSize_) : 0;
			if (ret != 0) {
				std::cerr<<__FILE__;
				std::cerr<<"\tpthread_attr_setstacksize() has returned non-zero "<<ret<<"\n";
				std::cerr<<"\tIt is possible (but no certain) that the following error";
				std::cerr<<"\thappened.\n";
				std::cerr<<"\tEINVAL The stack size is less than ";
				std::cerr<<"PTHREAD_STACK_MIN (16384) bytes.\n";
				std::cerr<<"\tI will ignore this error and let you continue\n";
			}

			ret = pthread_attr_init(attr[j]);
			checkForError(ret);

			ret = pthread_create(&thread_id[j],
			                     attr[j],
			                     mpiFunctionWrapper<SomeLambdaType, LoadBalancerType>,
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
	size_t stackSize_;
};
}
#endif // INTER_NODE_MPI_H
