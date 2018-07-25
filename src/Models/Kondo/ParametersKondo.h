#ifndef PARAMETERSKONDO_H
#define PARAMETERSKONDO_H

#include "ParametersModelBase.h"

namespace Dmrg {
template<typename RealType, typename QnType>
struct ParametersKondo : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	template<typename IoInputType>
	ParametersKondo(IoInputType& io) : BaseType(io, false)
	{
		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		io.readline(twiceTheSpin,"HeisenbergTwiceS=");
		if (twiceTheSpin != 1)
			err("ParametersKondo accepts only HeisenbergTwiceS=1 for now\n");
		io.read(potentialV, "potentialV");
		io.read(hubbardU, "hubbardU");
		io.read(kondoJ, "kondoJ");
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersKondo";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/kondoJ", kondoJ);
	}

	SizeType twiceTheSpin;
	VectorRealType potentialV;
	VectorRealType hubbardU;
	VectorRealType kondoJ;
};
}
#endif // PARAMETERSKONDO_H
