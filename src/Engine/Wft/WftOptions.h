#ifndef WFTOPTIONS_H
#define WFTOPTIONS_H
#include "Complex.h"

namespace Dmrg {

template<typename VectorWithOffsetType_>
struct WftOptions {

	typedef typename VectorWithOffsetType_::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	enum AccelEnum {ACCEL_NONE, ACCEL_PATCHES, ACCEL_BLOCKS, ACCEL_SVD};

	WftOptions(ProgramGlobals::DirectionEnum dir1,
	           PsimagLite::String options,
	           bool f,
	           bool b,
	           RealType d)
	    : twoSiteDmrg(options.find("twositedmrg") != PsimagLite::String::npos),
	      kronLoadBalance(options.find("KronLoadBalance") != PsimagLite::String::npos),
	      firstCall(f),
	      bounce(b),
	      dir(dir1),
	      accel((twoSiteDmrg) ? ACCEL_BLOCKS : ACCEL_PATCHES),
	      denseSparseThreshold(d)
	{
		if (options.find("wftAccelSvd") != PsimagLite::String::npos)
			accel = ACCEL_SVD;

		if (options.find("wftNoAccel") != PsimagLite::String::npos)
			accel = ACCEL_NONE;

		if (accel == ACCEL_SVD && twoSiteDmrg)
			err("wftAccelSvd not yet supported with twositedmrg\n");
	}

	void read(PsimagLite::IoSelector::In& io, PsimagLite::String label)
	{
		io.read(dir, label + "/dir");
		io.read(twoSiteDmrg, label + "/twoSiteDmrg");
		io.read(accel, label + "/accel");
		io.read(kronLoadBalance, label + "/kronLoadBalance");
		io.read(firstCall, label + "/firstCall");
		io.read(bounce, label + "/bounce");
		io.read(denseSparseThreshold, label + "/denseSparseThreshold");
	}

	void write(PsimagLite::IoSelector::Out& io, PsimagLite::String label) const
	{
		io.createGroup(label);
		io.write(dir, label + "/dir");
		io.write(twoSiteDmrg, label + "/twoSiteDmrg");
		io.write(accel, label + "/accel");
		io.write(kronLoadBalance, label + "/kronLoadBalance");
		io.write(firstCall, label + "/firstCall");
		io.write(bounce, label + "/bounce");
		io.write(denseSparseThreshold, label + "/denseSparseThreshold");
	}

	bool twoSiteDmrg;
	bool kronLoadBalance;
	bool firstCall;
	bool bounce;
	ProgramGlobals::DirectionEnum dir;
	AccelEnum accel;
	RealType denseSparseThreshold;

private:

	void accelMustBeNone(SizeType x) const
	{
		if (accel != ACCEL_NONE) {
			err("WFTOptions: More than one accel mode specified. Specify none or 1\n");
			return;
		}

		if (x == 0 || twoSiteDmrg) return;

		err("WFTOptions: onesitedmrg only with ACCEL_NONE or ACCEL_PATCHES\n");
	}
};

}
#endif // WFTOPTIONS_H
