#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"

namespace Dmrg {

template<typename ModelHelperBase>
class LinkProductBase {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	template<typename SomeStructType>
	static void connectorDofs(VectorSizeType& edofs,
	                          SizeType,
	                          SizeType,
	                          const SomeStructType&)
	{
		edofs[0] = edofs[1] = 0;
	}

	static SizeType dofsAllocationSize() { return 2; }

};
}
#endif // LINKPRODUCTBASE_H
