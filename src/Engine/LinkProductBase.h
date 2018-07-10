#ifndef LINKPRODUCTBASE_H
#define LINKPRODUCTBASE_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename ModelHelperType, typename GeometryType>
class LinkProductBase {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename ModelHelperType::RealType RealType;

	LinkProductBases(SizeType orbitals, bool hot, bool hasSpinOrbit)
	    : orbitals_(orbitals), hot_(hot), hasSpinOrbit_(hasSpinOrbit)
	{}

	bool hasSpinOrbit() const { return hasSpinOrbit_; }

	bool hot() const { return hot_; }

	const SizeType& orbitals() const { return orbitals_; }

	virtual SizeType dofs(SizeType, const AdditionalDataType&) = 0;

	virtual SizeType terms()  = 0;

	virtual void setLinkData(SizeType,
	                         SizeType,
	                         bool,
	                         ProgramGlobals::FermionOrBosonEnum&,
	                         PairType&,
	                         std::pair<char,char>&,
	                         SizeType&,
	                         RealType&,
	                         SizeType&,
	                         const AdditionalDataType&)  = 0;

	virtual void connectorDofs(VectorSizeType& edofs,
	                           SizeType,
	                           SizeType,
	                           const AdditionalDataType&)
	{
		edofs[0] = edofs[1] = 0;
	}

	virtual SizeType dofsAllocationSize() { return 2; }

	virtual void valueModifier(ComplexOrRealType&,
	                           SizeType,
	                           SizeType,
	                           bool,
	                           const AdditionalDataType&)
	{}

private:

	SizeType orbitals_;
	bool hot_;
	bool hasSpinOrbit_;
};
}
#endif // LINKPRODUCTBASE_H
