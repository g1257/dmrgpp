#ifndef PARAMETERSMODELBASE_H
#define PARAMETERSMODELBASE_H
#include "Io/IoNg.h"

namespace Dmrg {

template<typename RealType, typename QnType>
class ParametersModelBase {

public:

	template<typename IoInputType>
	ParametersModelBase(IoInputType&, bool)
	{}

	void write(PsimagLite::String,
	           PsimagLite::IoNg::Out::Serializer&) const {}
};
}
#endif // PARAMETERSMODELBASE_H
