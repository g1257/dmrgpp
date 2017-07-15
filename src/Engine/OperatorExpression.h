#ifndef OPERATOREXPRESSION_H
#define OPERATOREXPRESSION_H
#include "OperatorSpec.h"

namespace Dmrg {

template<typename ModelType>
class OperatorExpression {

	typedef OperatorSpec<ModelType> OperatorSpecType;
	typedef typename ModelType::OperatorType OperatorType;

public:

	OperatorExpression(const ModelType& model)
	    : operatorSpec_(model)
	{}

	// make sure that if a site is specified in an opsec, it
	// is the same in all others
	OperatorType operator()(PsimagLite::String opLabel, int& site2)
	{
		return operatorSpec_(opLabel, site2);
	}

private:

	OperatorSpecType operatorSpec_;
};
}
#endif // OPERATOREXPRESSION_H
