#ifndef LOADBALANCERDEFAULT_H
#define LOADBALANCERDEFAULT_H
#include "Vector.h"

namespace PsimagLite {

class LoadBalancerDefault {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	LoadBalancerDefault(const VectorSizeType& weights, SizeType nthreads)
	    : taskNumber_(nthreads)
	{
		VectorSizeType workLoad(nthreads, 0);

		SizeType ntasks = weights.size();
		for (SizeType i = 0; i < ntasks; ++i) {
			SizeType thread = findThreadWithLightestWork(workLoad);
			// assign work to thread
			assert(thread < taskNumber_.size());
			taskNumber_[thread].push_back(i);
			// update work loads
			assert(thread < workLoad.size());
			workLoad[thread] += weights[i];
		}

#ifndef NDEBUG
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
}; // class LoadBalancerDefault
} // namespace PsimagLite
#endif // LOADBALANCERDEFAULT_H
