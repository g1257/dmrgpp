#ifndef LOADBALANCER_MPI_H
#define LOADBALANCER_MPI_H
#include "Sort.h"
#include "Vector.h"

namespace PsimagLite {

class LoadBalancerMpi {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	LoadBalancerMpi(SizeType ntasks, SizeType nthreads)
	    : blockSize_(ntasks / nthreads)
	{
		if (ntasks < nthreads && ntasks > 0) {
			nthreads = ntasks;
			blockSize_ = 1;
		}

		assert(nthreads > 0);
		if ((ntasks % nthreads) != 0)
			++blockSize_;
	}

	SizeType blockSize(SizeType) const { return blockSize_; }

	SizeType taskNumber(SizeType threadNum, SizeType p) const
	{
		return p + threadNum * blockSize_;
	}

private:

	SizeType blockSize_;
}; // class LoadBalancerMpi
} // namespace PsimagLite
#endif // LOADBALANCER_MPI_H
