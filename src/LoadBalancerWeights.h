#ifndef LOADBALANCERWEIGHTS_H
#define LOADBALANCERWEIGHTS_H
#include "Vector.h"
#include "Sort.h"

namespace PsimagLite {

class LoadBalancerWeights {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	LoadBalancerWeights(SizeType ntasks, SizeType nthreads)
	    : LoadBalancerWeights(VectorSizeType(ntasks, 1), nthreads) // ctor delegation
	{}

	LoadBalancerWeights(const VectorSizeType& weights, SizeType nthreads)
	    : taskNumber_(nthreads)
	{
		SizeType ntasks = weights.size();
		if (ntasks < nthreads && ntasks > 0) nthreads = ntasks;
		VectorSizeType workLoad(nthreads, 0);
		VectorSizeType weights2 = weights;
		VectorSizeType iperm(ntasks, 0);
		Sort<VectorSizeType> sort;
		sort.sort(weights2,iperm);

		for (SizeType iii = 0; iii < ntasks; ++iii) {
			SizeType ii = ntasks - 1 - iii; // because sort is ascending
			SizeType thread = findThreadWithLightestWork(workLoad);
			// assign work to thread
			assert(thread < taskNumber_.size());
			taskNumber_[thread].push_back(iperm[ii]);
			// update work loads
			assert(thread < workLoad.size()); 
			workLoad[thread] += weights[iperm[ii]];
		}

#ifdef DEBUG_PTHREADS_NG
		for (SizeType i = 0; i < nthreads; ++i) {
			SizeType n = taskNumber_[i].size();
			std::cout<<n<<" Indices allocated to thread "<<i<<": ";
			for (SizeType j = 0; j < n; ++j)
				std::cout<<taskNumber_[i][j]<<" ";
			std::cout<<"\n";
		}
#endif
	}

	SizeType blockSize(SizeType threadNum) const
	{
		assert(threadNum < taskNumber_.size());
		return taskNumber_[threadNum].size();
	}

	int taskNumber(SizeType threadNum, SizeType p) const
	{
		assert(threadNum < taskNumber_.size());
		assert(p < taskNumber_[threadNum].size());
		return taskNumber_[threadNum][p];
	}

private:

	SizeType findThreadWithLightestWork(const VectorSizeType& workLoad) const
	{
		return std::min_element(workLoad.begin(), workLoad.end()) - workLoad.begin();
	}

	PsimagLite::Vector<VectorSizeType>::Type taskNumber_;
}; // class LoadBalancerWeights
} // namespace PsimagLite
#endif // LOADBALANCERWEIGHTS_H
