#ifndef DUMPHAMILTONIAN_HH
#define DUMPHAMILTONIAN_HH
#include "CrsMatrix.h"
#include "LoopSiteDirection.hh"
#include "PredicateAwesome.h"
#include "ProgramGlobals.h"
#include "Vector.h"
#include <fstream>
#include <string>

namespace Dmrg
{

template <typename ParametersDmrgType, typename ComplexOrRealType>
class DumpHamiltonian
{

public:

	using VectorType = std::vector<ComplexOrRealType>;
	using SparseMatrixType = PsimagLite::CrsMatrix<ComplexOrRealType>;

	DumpHamiltonian(const ParametersDmrgType& parameters)
	    : parameters_(parameters)
	    , predicate_(parameters.dumpHamiltonian)
	    , enabled_(parameters.options.isSet("debugmatrix") && !predicate_.empty())
	    , counter_(0)
	    , filename_(ProgramGlobals::rootName(parameters.filename) + "_dumpH.txt")
	{
		if (!enabled_) {
			return;
		}

		fout_.open(filename_.c_str());
		fout_.precision(parameters_.precision);
		fout_ << "DumpHamiltonian for DMRG++ version " << DMRGPP_VERSION << "\n";
		fout_ << "All matrices below are in MM Format\n";
		fout_ << buildMMFormatInfo();
	}

	~DumpHamiltonian()
	{
		if (!enabled_) {
			return;
		}

		fout_.close();
		std::cout << "DumpHamiltonian: file " << filename_ << " has been written.\n";
	}

	void save(const SparseMatrixType& crs, const VectorType& initial_vector, const LoopSiteDirection& loop_site_dir)
	{
		if (!enabled_) {
			return;
		}

		if (loop_site_dir.direction == ProgramGlobals::DirectionEnum::INFINITE) {
			return;
		}

		std::string predicate = predicate_;
		PsimagLite::replaceAll(predicate, "s", ttos(loop_site_dir.site));
		PsimagLite::PredicateAwesome<> pAwesome(predicate);
		bool b = pAwesome.isTrue("l", loop_site_dir.loopIndex);

		if (!b) {
			return;
		}

		fout_ << "Instance=" << counter_ << "\n";
		fout_ << "Loop=" << loop_site_dir.loopIndex << "\n";
		fout_ << "Site=" << loop_site_dir.site << "\n";
		fout_ << "Matrix_in_MM_Format_follows\n";

		printInMMFormat(crs);

		fout_ << "InitVector_follows\n";
		fout_ << initial_vector;
		++counter_;
	}

private:

	void printInMMFormat(const SparseMatrixType& crs)
	{
		SizeType rows = crs.rows();
		SizeType cols = crs.cols();
		SizeType nonzeroes = crs.nonZeros();
		fout_ << buildMMFormatHeader() << "\n";
		fout_ << rows << " " << cols << " " << nonzeroes << "\n";
		for (SizeType i = 0; i < rows; ++i) {
			SizeType start = crs.getRowPtr(i);
			SizeType end = crs.getRowPtr(i + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType col = crs.getCol(k);
				const ComplexOrRealType& value = crs.getValue(k);
				// one-based output required by MM Format
				fout_ << (i + 1) << " " << (col + 1) << " " << value << "\n";
			}
		}
	}

	std::string buildMMFormatInfo()
	{
		return std::string("%=================================================================================\n"
				   "%\n"
				   "% A sparse MxN matrix with L \n"
				   "% nonzeros in the Matrix Market format is represented as follows:\n"
				   "%\n"
				   "% +----------------------------------------------+\n"
				   "% |%%MatrixMarket matrix coordinate real general | <--- header line\n"
				   "% |%                                             | <--+\n"
				   "% |% comments                                    |    |-- 0 or more comment lines\n"
				   "% |%                                             | <--+        \n"
				   "% |    M  N  L                                   | <--- rows, columns, entries\n"
				   "% |    I1  J1  A(I1, J1)                         | <--+\n"
				   "% |    I2  J2  A(I2, J2)                         |    |\n"
				   "% |    I3  J3  A(I3, J3)                         |    |-- L lines\n"
				   "% |        . . .                                 |    |\n"
				   "% |    IL JL  A(IL, JL)                          | <--+\n"
				   "% +----------------------------------------------+   \n"
				   "%\n"
				   "% Indices are 1-based, i.e. A(1,1) is the first element.\n"
				   "%\n"
				   "%=================================================================================\n");
	}

	std::string buildMMFormatHeader()
	{
		return std::string("%%MatrixMarket matrix coordinate real general");
	}

	const ParametersDmrgType& parameters_;
	std::string predicate_;
	SizeType counter_;
	std::string filename_;
	bool enabled_;
	std::ofstream fout_;
};
}
#endif // DUMPHAMILTONIAN_HH
