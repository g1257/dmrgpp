#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelHelperType_, typename GeometryType_>
class LinkProductBase {

public:

	typedef ModelHelperType_ ModelHelperType;
	typedef GeometryType_ GeometryType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename ModelHelperType::RealType RealType;
	typedef std::pair<SizeType, SizeType> PairSizeType;

	virtual ~LinkProductBase() {}

	virtual SizeType terms() const = 0;

	virtual SizeType dofs(SizeType, const AdditionalDataType&) const = 0;

	virtual void setLinkData(SizeType,
	                         SizeType,
	                         bool,
	                         ProgramGlobals::FermionOrBosonEnum&,
	                         PairSizeType&,
	                         std::pair<char,char>&,
	                         SizeType&,
	                         RealType&,
	                         SizeType&,
	                         const AdditionalDataType&) const = 0;

	virtual void connectorDofs(VectorSizeType& edofs,
	                           SizeType,
	                           SizeType,
	                           const AdditionalDataType&) const
	{
		edofs[0] = edofs[1] = 0;
	}

	virtual SizeType dofsAllocationSize() const { return 2; }

	virtual void valueModifier(ComplexOrRealType&,
	                           SizeType,
	                           SizeType,
	                           bool,
	                           const AdditionalDataType&) const
	{}
};
}
#endif // LINKPRODUCTBASE_H
