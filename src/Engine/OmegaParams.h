#ifndef OMEGAPARAMS_H
#define OMEGAPARAMS_H
#include "Vector.h"
#include "InputCheck.h"

namespace Dmrg {

template<typename InputNgType, typename RealType_>
class OmegaParams {

public:

	typedef RealType_ RealType;

	OmegaParams(PsimagLite::String data)
	{
		Dmrg::InputCheck inputCheck;
		typename InputNgType::Writeable ioWriteable(inputCheck, data);
		typename InputNgType::Readable io(ioWriteable);
		configure(io);
	}

	OmegaParams(typename InputNgType::Readable& io)
	{
		configure(io);
	}

	void configure(typename InputNgType::Readable& io)
	{
		io.readline(begin_, "OmegaBegin=");
		io.readline(step_, "OmegaStep=");
		io.readline(total_, "OmegaTotal=");
		io.readline(offset_, "OmegaOffset=");
		io.readline(obs_, "Observable=");
	}

	RealType omega(SizeType i) const
	{
		return i*step_ + begin_;
	}

	PsimagLite::String observable() const { return obs_; }

	SizeType offset() const { return offset_; }

	SizeType total() const { return total_; }

private:

	RealType begin_;
	RealType step_;
	SizeType offset_;
	SizeType total_;
	PsimagLite::String obs_;
};

}
#endif // OMEGAPARAMS_H
