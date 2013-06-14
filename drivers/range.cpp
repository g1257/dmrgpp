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
	    : sum_(nthreads),v_(total)
	{}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      typename ConcurrencyType::MutexType* myMutex)
	{
		for (SizeType p=0;p<blockSize;p++) {
			SizeType taskNumber = threadNum*blockSize + p;
			if (taskNumber>=total) break;

			std::cout<<"This is thread number "<<threadNum;
			std::cout<<" and taskNumber="<<taskNumber<<"\n";

			SizeType ind = ConcurrencyType::storageIndex(threadNum);
			sum_[ind] += taskNumber;
			v_[taskNumber] = taskNumber * taskNumber;
		}
	}

	template<typename SomeParallelType>
	SizeType sum(SomeParallelType& p)
	{
		if (ConcurrencyType::mode == ConcurrencyType::MPI) {
			SizeType tmp = sum_[0];
			p.gather(tmp);
			p.bcast(tmp);
			return tmp;
		}
		return sumInternal(sum_);
	}

	template<typename SomeParallelType>
	const PsimagLite::Vector<SizeType>::Type& v(SomeParallelType& p)
	{
		if (ConcurrencyType::mode == ConcurrencyType::MPI) {
			p.allGather(v_);
		}
		return v_;
	}

private:

	template<typename SomeVectorType>
	typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
	                              typename SomeVectorType::value_type>::Type
	sumInternal(SomeVectorType& v)
	{
		typename SomeVectorType::value_type tmp = 0;
		for (size_t i=0;i<v.size();i++) {
			tmp += v[i];
		}
		return tmp;
	}

	PsimagLite::Vector<SizeType>::Type sum_;
	PsimagLite::Vector<SizeType>::Type v_;
};

int main(int argc,char *argv[])
{

	typedef PsimagLite::Concurrency ConcurrencyType;
	ConcurrencyType concurrency(argc,argv);

	typedef MyLoop HelperType;
	typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
	ParallelizerType threadObject;

	if (argc < 3) return 1;

	SizeType total = atoi(argv[1]);

	SizeType nthreads = 1;
	if (ConcurrencyType::mode == ConcurrencyType::PTHREADS) {
		if (argc != 3) return 1;
		nthreads = atoi(argv[2]);
	}
	ParallelizerType::setThreads(nthreads);

	HelperType helper(nthreads,total);

	std::cout<<"Using "<<threadObject.name();
	std::cout<<" with "<<threadObject.threads()<<" threads.\n";
	threadObject.loopCreate(total,helper);

	SizeType sum = helper.sum(threadObject);

	if (concurrency.root()) std::cout<<"sum="<<sum<<"\n";

	std::cout<<helper.v(threadObject);
}
