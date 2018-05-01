/** How to compile and run this driver
 *
 * Serial version:
 *
 * g++ -g3 -DNDEBUG  -Werror -Wall -I../src  \
 *  -lm  -lpthread   range.cpp   -o range \
 * -lblas -llapack
 *
 * And run it with:
 *
 * ./range
 *
 * Parallel version:
 *
 * mpicxx -DUSE_MPI -g3 -DNDEBUG  -Werror -Wall -I../src \
 *  -lm  -lpthread   range.cpp   -o range \
 * -lblas -llapack
 *
 * And run it with:
 *
 * your batch system script
 *
 */
#define USE_PTHREADS_OR_NOT_NG
#include "Concurrency.h"
#include "Parallelizer.h"
#include <iostream>
#include <unistd.h>

class MyLoop {

	typedef PsimagLite::Concurrency ConcurrencyType;

public:

	MyLoop(SizeType nthreads, SizeType total)
	    : sum_(ConcurrencyType::storageSize(nthreads)),
	      v_(total,0)
	{}

	SizeType tasks() const { return v_.size(); }

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		sleep(1);

		SizeType ind = ConcurrencyType::storageIndex(threadNum);
		sum_[ind] += taskNumber;
		v_[taskNumber] = taskNumber * taskNumber;
	}

	void sync()
	{
		if (ConcurrencyType::hasPthreads()) {
			SizeType tmp = PsimagLite::sum(sum_);
			sum_[0] = tmp;
		}

		if (ConcurrencyType::hasMpi()) {
			SizeType tmp = sum_[0];
			PsimagLite::MPI::allReduce(tmp);
			sum_[0] = tmp;

			PsimagLite::MPI::allReduce(v_);
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

	ConcurrencyType concurrency(&argc,&argv,nthreads);

	typedef MyLoop HelperType;
	typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
	ParallelizerType threadObject(ConcurrencyType::codeSectionParams);

	if (argc < 2) return 1;

	SizeType total = atoi(argv[1]);

	HelperType helper(nthreads,total);

	threadObject.loopCreate(helper);

	helper.sync();

	SizeType sum = helper.sum();

	if (ConcurrencyType::root()) {
		std::cout<<"Using "<<threadObject.name()<<" mode= "<<ConcurrencyType::mode;
		std::cout<<" with "<<threadObject.threads();
		std::cout<<" threads or mpi procs.\n";
		std::cout<<"sum="<<sum<<"\n";
		std::cout<<helper.v();
	}

}
