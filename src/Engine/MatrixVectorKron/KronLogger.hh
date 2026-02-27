#ifndef KRONLOGGER_HH
#define KRONLOGGER_HH
#include "Basis.h"
#include "BasisWithOperators.h"
#include "CrsMatrix.h"
#include "InitKronHamiltonian.h"
#include "LeftRightSuper.h"
#include "MatrixMarket.hh"
#include "ProgressIndicator.h"
#include "PsimagLite.h"
#include <fstream>

namespace Dmrg {

template <typename ModelType> class KronLogger {
public:

	using InitKronType            = InitKronHamiltonian<ModelType>;
	using ArrayOfMatStructType    = typename InitKronType::ArrayOfMatStructType;
	using MatrixDenseOrSparseType = typename ArrayOfMatStructType::MatrixDenseOrSparseType;
	using ComplexOrRealType       = typename InitKronType::ComplexOrRealType;
	using MatrixMarketType        = MatrixMarket<ComplexOrRealType>;

	KronLogger(const InitKronType& init_kron)
	    : progress_("KronLogger")
	    , fout_(nullptr)
	    , init_kron_(init_kron)
	{
		if (!init_kron.params().options.isSet("KroneckerDumper")) {
			return;
		}

		if (counter_ >= init_kron.params().dumperEnd) {
			return;
		}

		if (counter_++ < init_kron.params().dumperBegin) {
			return;
		}

		std::string filename = buildFilename(init_kron.params().filename, counter_);
		fout_                = new std::ofstream(filename);
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

		SizeType nC      = init_kron_.connections();
		SizeType total   = init_kron_.numberOfPatches(InitKronType::OLD);
		SizeType offsetX = init_kron_.offsetForPatches(InitKronType::NEW, outPatch);
		*fout_ << "outPatch=" << outPatch << "\n";
		*fout_ << "nC=" << nC << " total=" << total << " offsetX=" << offsetX << "\n";
	}

	void two(SizeType inPatch)
	{
		if (!fout_)
			return;

		SizeType offsetY = init_kron_.offsetForPatches(InitKronType::OLD, inPatch);
		*fout_ << separationLevel(1) << "inPatch=" << inPatch << "\n";
		*fout_ << separationLevel(1) << "offsetY=" << offsetY << "\n";
	}

	void three(SizeType outPatch, SizeType inPatch, SizeType ic)
	{
		if (!fout_)
			return;

		const ArrayOfMatStructType& xiStruct = init_kron_.xc(ic);
		const ArrayOfMatStructType& yiStruct = init_kron_.yc(ic);

		const bool performTranspose = (init_kron_.useLowerPart() && (outPatch < inPatch));

		const MatrixDenseOrSparseType* Amat
		    = performTranspose ? xiStruct(inPatch, outPatch) : xiStruct(outPatch, inPatch);

		const MatrixDenseOrSparseType* Bmat
		    = performTranspose ? yiStruct(inPatch, outPatch) : yiStruct(outPatch, inPatch);

		*fout_ << separationLevel(2) << "ic=" << ic << "\n";
		*fout_ << separationLevel(2) << "performTranspose=" << performTranspose << "\n";

		// TODO: Convert printing to Matrix Market Format
		assert(Amat);
		*fout_ << separationLevel(2)
		       << "Matrix A follows in format: " + matrixFormat(Amat->isDense()) + "\n";
		printMatrixDenseOrSparse(*Amat);

		assert(Bmat);
		*fout_ << separationLevel(2)
		       << "Matrix B follows in format: " + matrixFormat(Bmat->isDense()) + "\n";
		printMatrixDenseOrSparse(*Bmat);
	}

private:

	static std::string separationLevel(SizeType n)
	{
		std::string tmp;
		for (SizeType i = 0; i < n; ++i) {
			tmp += " ";
		}

		return tmp;
	}

	static std::string matrixFormat(bool b) { return (b) ? "dense" : "sparse"; }

	void printMatrixDenseOrSparse(const MatrixDenseOrSparseType& mat)
	{
		if (!fout_) {
			err("printMatrixDenseOrSparse: InternalError: KronLogger not enabled!?\n");
		}

		if (mat.isDense()) {
			*fout_ << mat.dense();
		} else {
			MatrixMarketType matrix_market(mat.sparse());
			matrix_market.print(*fout_);
		}
	}

	static std::string buildFilename(const std::string& filename, SizeType n)
	{
		// FIXME find the last . in filename if any
		size_t dot_index = filename.find(".");

		std::string root = filename.substr(0, dot_index);
		return "kron_" + root + ttos(n) + ".txt";
	}

	static SizeType               counter_;
	PsimagLite::ProgressIndicator progress_;
	std::ofstream*                fout_;
	const InitKronType&           init_kron_;
};

template <typename ModelType> SizeType KronLogger<ModelType>::counter_ = 0;
}
#endif // KRONLOGGER_HH
