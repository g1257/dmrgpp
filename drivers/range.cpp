/** How to compile and run this driver
 * 
 * Serial version:
 * 
 * g++ -g3 -DNDEBUG  -Werror -Wall -I../src -I../src/JSON \
 * -I../src/JSON/JsonParser -lm  -lpthread   range.cpp   -o range \
 * /usr/lib64/libblas.so.3 /usr/lib64/liblapack.so.3
 * 
 * And run it with:
 * 
 * ./range
 * 
 * Parallel version:
 * 
 * mpicxx -DUSE_MPI -g3 -DNDEBUG  -Werror -Wall -I../src -I../src/JSON \
 * -I../src/JSON/JsonParser -lm  -lpthread   range.cpp   -o range \
 * /usr/lib64/libblas.so.3 /usr/lib64/liblapack.so.3
 * 
 * And run it with:
 * 
 * your batch system script
 *
 */

#include "Concurrency.h"
#include "Parallelizer.h"
#include <iostream>

class MyLoop {

	typedef PsimagLite::Concurrency ConcurrencyType;

public:

	MyLoop(SizeType nthreads,SizeType total)
	    : sum_(ConcurrencyType::storageSize(nthreads)),
	      v_(total,0)
	{}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      ConcurrencyType::MutexType* myMutex)
	{
		for (SizeType p=0;p<blockSize;p++) {
			SizeType taskNumber = threadNum*blockSize + p;
			if (taskNumber>=total) break;

//			std::cout<<"This is thread number "<<threadNum;
//			std::cout<<" and taskNumber="<<taskNumber<<"\n";

			SizeType ind = ConcurrencyType::storageIndex(threadNum);
			sum_[ind] += taskNumber;
			v_[taskNumber] = taskNumber * taskNumber;
		}
	}

	template<typename SomeParallelType>
	void sync(SomeParallelType& p)
	{
		if (ConcurrencyType::mode == ConcurrencyType::MPI) {
			assert(sum_.size() == 1);
			SizeType tmp = sum_[0];
			p.allReduce(tmp);
			sum_[0] = tmp;

			p.allReduce(v_);
		} else {
			SizeType tmp = PsimagLite::sum(sum_);
			sum_[0] = tmp;
		}
	}

	SizeType sum() const
	{
		assert(sum_.size() > 0);
		return sum_[0];
	}

	const PsimagLite::Vector<SizeType>::Type& v() const
	{
		return v_;
	}

private:

	PsimagLite::Vector<SizeType>::Type sum_;
	PsimagLite::Vector<SizeType>::Type v_;
};

int main(int argc,char *argv[])
{

	typedef PsimagLite::Concurrency ConcurrencyType;

	SizeType nthreads = 1;
	if (argc == 3) nthreads = atoi(argv[2]);

	ConcurrencyType concurrency(argc,argv,nthreads);

	typedef MyLoop HelperType;
	typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);

	if (argc < 2) return 1;

	SizeType total = atoi(argv[1]);

	HelperType helper(nthreads,total);

	threadObject.loopCreate(total,helper);

	helper.sync(threadObject);

	SizeType sum = helper.sum();

	if (ConcurrencyType::root()) {
		std::cout<<"Using "<<threadObject.name();
		std::cout<<" with "<<threadObject.threads()<<" threads.\n";
		std::cout<<"sum="<<sum<<"\n";
		std::cout<<helper.v();
	}
		
}
