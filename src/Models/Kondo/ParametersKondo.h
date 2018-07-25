#ifndef PARAMETERSKONDO_H
#define PARAMETERSKONDO_H

#include "ParametersModelBase.h"

namespace Dmrg {
template<typename RealType, typename QnType>
struct ParametersKondo : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;

	template<typename IoInputType>
	ParametersKondo(IoInputType& io) : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		io.readline(twiceTheSpin,"HeisenbergTwiceS=");
		if (twiceTheSpin != 1)
			err("ParametersKondo accepts only HeisenbergTwiceS=1 for now\n");
	}

	SizeType twiceTheSpin;
};
}
#endif // PARAMETERSKONDO_H
