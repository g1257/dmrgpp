#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"

namespace Dmrg {

template<typename ModelHelperType, typename GeometryType>
class LinkProductBase {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;

	virtual void connectorDofs(VectorSizeType& edofs,
	                          SizeType,
	                          SizeType,
	                          const AdditionalDataType&)
	{
		edofs[0] = edofs[1] = 0;
	}

	virtual SizeType dofsAllocationSize() { return 2; }

	virtual void valueModifier(SparseElementType&,
	                          SizeType,
	                          SizeType,
	                          bool,
	                          const AdditionalDataType&)
	{}
};
}
#endif // LINKPRODUCTBASE_H
