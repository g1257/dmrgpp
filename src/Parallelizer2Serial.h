#ifndef PARALLELIZER2SERIAL_H
#define PARALLELIZER2SERIAL_H
#include "Vector.h"
#include "CodeSectionParams.h"

namespace PsimagLite {

template<typename = int>
class Parallelizer2 {

public:

	Parallelizer2(const CodeSectionParams& codeParams)
	{
		if (codeParams.npthreads != 1)
			throw RuntimeError("Please compile with -DUSE_PTHREADS\n");
	}

	template<typename SomeLambdaType>
	void parallelFor(SizeType start, SizeType end, const SomeLambdaType& lambda)
	{
		for (SizeType i = start; i < end; ++i)
			lambda(i, 0);
	}

	SizeType numberOfThreads() const { return 1; }

	String name() const { return "serial"; }
};
}
#endif // PARALLELIZER2SERIAL_H
