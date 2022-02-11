#ifndef TRUNCATIONCONTROL_H
#define TRUNCATIONCONTROL_H
#include "PsimagLite.h"
#include "ProgressIndicator.h"
#include "InputNg.h"
#include "InputCheck.h"

namespace Dmrg {

template<typename RealType>
class TruncationControl {

public:

	using VectorStringType = PsimagLite::Vector<PsimagLite::String>::Type;
	using InputNgReadableType = PsimagLite::InputNg<InputCheck>::Readable;

	TruncationControl() : tolerance_(-1), mMin_(0) {}

	void read(InputNgReadableType& io, SizeType keptStatesInfinite, bool hasTwoSiteDmrg)
	{
		mMin_ = keptStatesInfinite;

		try {
			PsimagLite::String s;
			VectorStringType tokens;
			io.readline(s, "TruncationTolerance=");
			PsimagLite::split(tokens, s, ",");
			tolerance_ = PsimagLite::atof(tokens[0].c_str());
			if (tokens.size() > 1)
				mMin_ = atoi(tokens[1].c_str());
			if (!hasTwoSiteDmrg) {
				std::cerr<<"WARNING: TruncationTolerance used without twositedmrg\n";
				std::cout<<"WARNING: TruncationTolerance used without twositedmrg\n";
			}
		} catch (std::exception&) {}
	}

	void write(PsimagLite::String label, PsimagLite::IoSerializer& ioSerializer) const
	{
		ioSerializer.createGroup(label);
		ioSerializer.write(label + "/tolerance", tolerance_);
		ioSerializer.write(label + "/mMin", mMin_);
	}

	void print(std::ostream& os, PsimagLite::ProgressIndicator& progress) const
	{
		if (tolerance_ < 0) return;
		PsimagLite::OstringStream msgg(os.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"has tolerance= "<<tolerance_;
		msg<<" minimum m= "<<mMin_;
		progress.printline(msgg, os);
	}

	void setTolerance(RealType tol) { tolerance_ = tol; }

	void setMmin(SizeType mMin) { mMin_ = mMin; }

	RealType tolerance() const { return tolerance_; }

	SizeType mMin() const { return mMin_; }

private:

	RealType tolerance_;
	SizeType mMin_;
};
}
#endif // TRUNCATIONCONTROL_H
