#ifndef INTER_NODE_MPI_H
#define INTER_NODE_MPI_H
#include "Mpi.h"
#include <iostream>
#include <algorithm>
#include "Vector.h"
#include <sched.h>
#include <unistd.h>
#include "TypeToString.h"
#include "LoadBalancerMpi.h"

namespace PsimagLite {

template<typename LoadBalancerType=LoadBalancerMpi>
class InterNode {

public:

	typedef LoadBalancerMpi::VectorSizeType VectorSizeType;

	InterNode(MPI::CommType comm)
	    : comm_(comm), mpiSize_(MPI::commSize(comm_)), mpiRank_(MPI::commRank(comm_))
	{}

	SizeType size() const { return mpiSize_; }

	String name() const { return "mpi"; }

	// no weights, no balancer ==> create weights, set all weigths to 1, delegate
	template<typename SomeLambdaType>
	void parallelFor(SizeType start, SizeType end, const SomeLambdaType& lambda)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(end - start, mpiSize_);
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
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights.size(), mpiSize_);
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
		SizeType blockSize = loadBalancer.blockSize(mpiRank_);

		for (SizeType p = 0; p < blockSize; ++p) {
			SizeType taskNumber = loadBalancer.taskNumber(mpiRank_, p);
			if (taskNumber + start >= end) break;
			lambda(taskNumber + start, mpiRank_);
		}

		MPI::barrier(comm_);
	}

private:

	MPI::CommType comm_;
	SizeType mpiSize_;
	SizeType mpiRank_;
};
}
#endif // INTER_NODE_MPI_H
