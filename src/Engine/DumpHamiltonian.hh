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
		fout_ << "CRS_Matrix_follows\n";
		fout_ << crs;
		fout_ << "InitVector_follows\n";
		fout_ << initial_vector;
		++counter_;
	}

private:

	const ParametersDmrgType& parameters_;
	std::string predicate_;
	SizeType counter_;
	std::string filename_;
	bool enabled_;
	std::ofstream fout_;
};
}
#endif // DUMPHAMILTONIAN_HH
