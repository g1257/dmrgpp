#ifndef DMRG_INTROSPECT_H
#define DMRG_INTROSPECT_H
#include "OperatorSpec.h"
#include "OptionsForIntrospect.hh"
#include <PsimagLite/CanonicalExpression.h>
#include <PsimagLite/PsimagLite.h>

namespace Dmrg {

class Introspect {
public:

	Introspect(bool enabled)
	    : enabled_(enabled)
	{ }

	template <typename ModelBaseType>
	bool operator()(const ModelBaseType& model, const OptionsForIntrospect& obsOptions)
	{
		if (!enabled_) {
			return false;
		}

		using ModelHelperType  = typename ModelBaseType::ModelHelperType;
		using OperatorsType    = typename ModelHelperType::OperatorsType;
		using OperatorType     = typename OperatorsType::OperatorType;
		using OperatorSpecType = Dmrg::OperatorSpec<ModelBaseType, OperatorType>;

		OperatorType                                      opC;
		const OperatorType                                opEmpty;
		OperatorSpecType                                  opSpec(model);
		int                                               site = -1;
		PsimagLite::CanonicalExpression<OperatorSpecType> canonicalExpression(opSpec);

		switch (obsOptions.introspect) {
		case OptionsForIntrospect::IntrospectEnum::EXPRESSION:
			canonicalExpression(opC, obsOptions.opexpr, opEmpty, site);
			opC.write(std::cout);
			break;
		case OptionsForIntrospect::IntrospectEnum::MODEL_BASIS:
			model.printBasis(obsOptions.site);
			break;
		case OptionsForIntrospect::IntrospectEnum::MODEL_HAMILTONIAN:
			model.printTerms();
			break;
		default:
			throw std::runtime_error("InternalError: operatorDriver\n");
		}

		return true;
	}

private:

	bool enabled_;
};
}
#endif // DMRG_INTROSPECT_H
