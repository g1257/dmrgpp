#ifndef PARALLELIZER2SERIAL_H
#define PARALLELIZER2SERIAL_H
#include "Vector.h"

namespace PsimagLite {

template<typename = int>
class Parallizer2 {

public:

	Parallizer2(SizeType threads)
	{
		if (threads != 1)
			throw RuntimeError("Please compile with -DUSE_PTHREADS\n");
	}

	template<typename SomeLambdaType>
	void parallelFor(const SomeLambdaType& lambda, SizeType n)
	{
		for (SizeType i = 0; i < n; ++i)
			lambda(i, 0);
	}
};
}
#endif // PARALLELIZER2SERIAL_H
