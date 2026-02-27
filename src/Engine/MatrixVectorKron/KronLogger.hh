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

	// using ComplexOrRealType = typename SuperBlockType::SparseElementType;
	// using SparseMatrixType = PsimagLite::CrsMatrix<ComplexOrRealType>;
	// using BasisType = Basis<SparseMatrixType>;
	// using BasisWithOperatorsType = BasisWithOperators<BasisType>;
	// using LeftRightSuperType = LeftRightSuper<BasisWithOperatorsType, SuperBlockType>;
	using InitKronType = InitKronHamiltonian<ModelType>;

	KronLogger(const InitKronType& init_kron, const std::string& filename)
	    : progress_("KronLogger")
	    , fout_(filename)
	    , init_kron_(init_kron)

	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "KronLogger: Hello from ctor\n";
		progress_.printline(msgg, fout_);
	}

	~KronLogger()
	{
		PsimagLite::OstringStream                     msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg << "KronLogger: Bye from dtor\n";
		progress_.printline(msgg, fout_);
	}

	void one(SizeType outPatch) { fout_ << "outPatch=" << outPatch << "\n"; }

	void two(SizeType inPatch) { fout_ << "->inPatch=" << inPatch << "\n"; }

	void three(SizeType ic) { fout_ << "-->ic=" << ic << "\n"; }

private:

	PsimagLite::ProgressIndicator progress_;
	std::ofstream                 fout_;
	const InitKronType&           init_kron_;
};
}
#endif // KRONLOGGER_HH
