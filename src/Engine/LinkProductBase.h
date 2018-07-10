#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"

namespace Dmrg {

template<typename ModelHelperType>
class LinkProductBase {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;

	template<typename SomeStructType>
	static void connectorDofs(VectorSizeType& edofs,
	                          SizeType,
	                          SizeType,
	                          const SomeStructType&)
	{
		edofs[0] = edofs[1] = 0;
	}

	static SizeType dofsAllocationSize() { return 2; }

	template<typename SomeStructType>
	static void valueModifier(SparseElementType&,
	                          SizeType,
	                          SizeType,
	                          bool,
	                          const SomeStructType&)
	{}


};
}
#endif // LINKPRODUCTBASE_H
