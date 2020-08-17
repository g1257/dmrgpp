#ifndef OMEGAPARAMS_H
#define OMEGAPARAMS_H
#include "Vector.h"
#include "InputCheck.h"

namespace Dmrg {

template<typename InputNgType, typename RealType>
struct OmegaParams {

	OmegaParams(PsimagLite::String data)
	{
		Dmrg::InputCheck inputCheck;
		typename InputNgType::Writeable ioWriteable(inputCheck, data);
		typename InputNgType::Readable io(ioWriteable);
		io.readline(begin, "OmegaBegin=");
		io.readline(step, "OmegaStep=");
		io.readline(total, "OmegaTotal=");
		io.readline(offset, "OmegaOffset=");
		io.readline(obs, "Observable=");
	}

	RealType begin;
	RealType step;
	SizeType offset;
	SizeType total;
	PsimagLite::String obs;
};

}
#endif // OMEGAPARAMS_H
