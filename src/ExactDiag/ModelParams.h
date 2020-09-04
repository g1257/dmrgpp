#ifndef CINC_MODELPARAMS_H
#define CINC_MODELPARAMS_H
#include "Vector.h"

namespace Dmft {

template<typename RealType>
struct ModelParams {

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	ModelParams(const VectorRealType& bathParams)
	{
		SizeType bath = bathParams.size()/2;
		assert((bathParams.size() & 1) == 0);
		sites = bath + 1;
	}

	SizeType sites;
};

}
#endif // CINC_MODELPARAMS_H
