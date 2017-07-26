#ifndef OPERATOREXPRESSION_H
#define OPERATOREXPRESSION_H
#include "OperatorSpec.h"
#include "PsimagLite.h"
#include "CanonicalExpression.h"

namespace Dmrg {

template<typename ModelType>
class OperatorExpression {

	typedef OperatorSpec<ModelType> OperatorSpecType;
	typedef typename ModelType::MyBasis BasisType;
	typedef typename ModelType::OperatorType OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::CanonicalExpression<OperatorSpecType> CanonicalExpressionType;

public:

	OperatorExpression(const ModelType& model)
	    : operatorSpec_(model), canonicalExpression_(operatorSpec_)
	{}

	// make sure that if a site is specified in an opsec, it
	// is the same in all others
	OperatorType operator()(PsimagLite::String opLabel, int& site2) const
	{
		return canonicalExpression_(opLabel, site2);
	}

private:

	OperatorSpecType operatorSpec_;
	CanonicalExpressionType canonicalExpression_;
};
}
#endif // OPERATOREXPRESSION_H
