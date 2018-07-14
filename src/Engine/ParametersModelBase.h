#ifndef PARAMETERSMODELBASE_H
#define PARAMETERSMODELBASE_H
#include "TargetQuantumElectrons.h"

namespace Dmrg {

template<typename RealType, typename QnType>
class ParametersModelBase {

public:

	template<typename IoInputType>
	ParametersModelBase(IoInputType& io, bool allowUpDown)
	    : targetQuantum_(io, allowUpDown)
	{}

	const TargetQuantumElectrons<RealType, QnType>& targetQuantum() const
	{
		return targetQuantum_;
	}

	void write(PsimagLite::String label,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		targetQuantum_.write(label, io);
	}

private:

	TargetQuantumElectrons<RealType, QnType> targetQuantum_;
};
}
#endif // PARAMETERSMODELBASE_H
