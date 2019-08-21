#ifndef LOADBALANCERDEFAULT_H
#define LOADBALANCERDEFAULT_H
#include "Vector.h"
#include "Sort.h"

namespace PsimagLite {

class LoadBalancerDefault {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	LoadBalancerDefault(SizeType ntasks, SizeType nthreads)
	    : blockSize_(ntasks/nthreads)
	{
		if ((ntasks % nthreads) != 0) ++blockSize_;
	}

	SizeType blockSize(SizeType) const
	{
		return blockSize_;
	}

	SizeType taskNumber(SizeType threadNum, SizeType p) const
	{
		return p + threadNum*blockSize_;
	}

private:

	SizeType blockSize_;
}; // class LoadBalancerDefault
} // namespace PsimagLite
#endif // LOADBALANCERDEFAULT_H
