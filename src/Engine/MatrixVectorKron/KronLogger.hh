#ifndef KRONLOGGER_HH
#define KRONLOGGER_HH
#include "Basis.h"
#include "BasisWithOperators.h"
#include "CrsMatrix.h"
#include "InitKronHamiltonian.h"
#include "LeftRightSuper.h"
#include "ProgressIndicator.h"
#include "PsimagLite.h"
#include <fstream>

namespace Dmrg {

template <typename ModelType> class KronLogger {
public:

	using InitKronType = InitKronHamiltonian<ModelType>;

	KronLogger(const InitKronType& init_kron, const std::string& filename)
	    : progress_("KronLogger")
	    , fout_(nullptr)
	    , init_kron_(init_kron)
	{
		if (!init_kron.params().options.isSet("kronlogger"))
			return;

		fout_ = new std::ofstream(filename);
		if (!fout_ or !fout_->good() or fout_->bad()) {
			delete fout_;
			fout_ = nullptr;
			err(std::string("Failed to create KronLogger file ") + filename + "\n");
		}

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "KronLogger: Hello from ctor\n";
		progress_.printline(msgg, *fout_);
	}

	~KronLogger()
	{
		if (!fout_)
			return;

		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "KronLogger: Bye from dtor\n";
		progress_.printline(msgg, *fout_);

		delete fout_;
		fout_ = nullptr;
	}

	void one(SizeType outPatch)
	{
		if (!fout_)
			return;

		*fout_ << "outPatch=" << outPatch << "\n";
	}

	void two(SizeType inPatch)
	{
		if (!fout_)
			return;

		*fout_ << "->inPatch=" << inPatch << "\n";
	}

	void three(SizeType ic)
	{
		if (!fout_)
			return;

		*fout_ << "-->ic=" << ic << "\n";
	}

private:

	PsimagLite::ProgressIndicator progress_;
	std::ofstream*                fout_;
	const InitKronType&           init_kron_;
};
}
#endif // KRONLOGGER_HH
