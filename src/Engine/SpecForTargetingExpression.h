#ifndef SPECFORTARGETINGEXPRESSION_H
#define SPECFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "OneOperatorSpec.h"
#include "CanonicalExpression.h"
#include <numeric>
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "PackIndices.h"
#include "AlgebraForTargetingExpression.h"

namespace Dmrg {

template<typename TargetingBaseType>
class SpecForTargetingExpression {

public:

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef AlgebraForTargetingExpression<TargetingBaseType> AlgebraType;
	typedef AlgebraType ResultType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename AlgebraType::AuxiliaryType AuxiliaryType;

	AlgebraType operator()(PsimagLite::String str, const AuxiliaryType& aux) const
	{
		return AlgebraType(str, aux);
	}

	static bool metaEqual(const ResultType&, const ResultType&) { return true; }

	static bool isEmpty(const ResultType& term)
	{
		return term.isEmpty();
	}
};
}
#endif // SPECFORTARGETINGEXPRESSION_H
