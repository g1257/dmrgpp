#ifndef PARAMETERSKONDO_H
#define PARAMETERSKONDO_H

#include "ParametersModelBase.h"

namespace Dmrg {
template <typename RealType, typename QnType>
struct ParametersKondo : public ParametersModelBase<RealType, QnType> {

	typedef ParametersModelBase<RealType, QnType> BaseType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	enum class ExtEnum
	{
		NONE,
		EXTENDED,
		EXTENDED2
	};

	template <typename IoInputType>
	ParametersKondo(IoInputType& io, PsimagLite::String option)
	    : BaseType(io, false)
	    , extended(ExtEnum::NONE)
	{
		if (option == "Ex") {
			extended = ExtEnum::EXTENDED;
			std::cout << "Ex detected\n";
		} else if (option == "Ex2") {
			extended = ExtEnum::EXTENDED2;
			std::cout << "Ex2 detected\n";
		} else if (option != "") {
			err("ParametersKondo: Invalid option " + option + "\n");
		}

		SizeType nsites = 0;
		io.readline(nsites, "TotalNumberOfSites=");
		io.readline(twiceTheSpin, "HeisenbergTwiceS=");
		if (twiceTheSpin != 1)
			err("ParametersKondo accepts only HeisenbergTwiceS=1 for now\n");
		io.read(potentialV, "potentialV");
		checkVector(potentialV, "potentialV", 2 * nsites);
		io.read(hubbardU, "hubbardU");
		checkVector(hubbardU, "hubbardU", nsites);
		io.read(kondoJ, "kondoJ");
		checkVector(kondoJ, "kondoJ", nsites);

		if (extended == ExtEnum::NONE)
			return;

		io.readline(kondoHx, "KondoHx=");
		io.readline(electronHx, "ElectronHx=");
		io.readline(pairingField, "PairingField=");
	}

	void write(PsimagLite::String label1, PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/ParametersKondo";
		io.createGroup(label);
		BaseType::write(label, io);
		io.write(label + "/twiceTheSpin", twiceTheSpin);
		io.write(label + "/potentialV", potentialV);
		io.write(label + "/hubbardU", hubbardU);
		io.write(label + "/kondoJ", kondoJ);
		io.write(label + "/extended", extended);

		if (extended == ExtEnum::NONE)
			return;

		io.write(label + "/kondoHx", kondoHx);
		io.write(label + "/electronHx", electronHx);
		io.write(label + "/pairingField", pairingField);
	}

	SizeType twiceTheSpin;
	VectorRealType potentialV;
	VectorRealType hubbardU;
	VectorRealType kondoJ;
	ExtEnum extended;
	RealType kondoHx;
	RealType electronHx;
	RealType pairingField;

private:

	void checkVector(const VectorRealType& v, PsimagLite::String str, SizeType n)
	{
		if (v.size() == n)
			return;
		err("Vector " + str + " should be of size " + ttos(n) + "\n");
	}
};
}
#endif // PARAMETERSKONDO_H
