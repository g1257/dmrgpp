#ifndef KRONECKERDUMPER_H
#define KRONECKERDUMPER_H
#include "Vector.h"
#include "TypeToString.h"
#include <fstream>
#include "../Version.h"

namespace Dmrg {

template<typename SparseMatrixType>
class KroneckerDumper {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	struct ParamsForKroneckerDumper {
		ParamsForKroneckerDumper(bool enable = false, SizeType b = 0, SizeType e = 0)
		    : enabled(enable),begin(b),end(e)
		{}

		bool enabled;
		SizeType begin;
		SizeType end;
	}; // struct ParamsForKroneckerDumper

	KroneckerDumper(const ParamsForKroneckerDumper* p)
	    : enabled_(p && p->enabled)
	{
		if (!enabled_) return;
		bool b = (p->end > 0 && counter_ >= p->end);
		if (counter_ < p->begin || b) {
			counter_++;
			enabled_ = false;
			return;
		}

		PsimagLite::String filename = "kroneckerDumper" + ttos(counter_) + ".txt";
		fout_.open(filename.c_str());
		fout_<<"#KroneckerDumper for DMRG++ version "<<DMRGPP_VERSION<<"\n";
	}

	~KroneckerDumper()
	{
		fout_<<"#EOF\n";
		fout_.close();
	}

	void push(const SparseMatrixType&,
	          const SparseMatrixType&,
	          const VectorBoolType&,
	          const VectorSizeType&,
	          SizeType start,
	          SizeType end) const
	{
		if (!enabled_) return;
	}

	void push(bool option,const SparseMatrixType& hamiltonian) const
	{
		if (!enabled_) return;
	}

private:

	static SizeType counter_;
	bool enabled_;
	std::ofstream fout_;
}; // class KroneckerDumpter

template<typename SparseMatrixType>
SizeType KroneckerDumper<SparseMatrixType>::counter_ = 0;

} // namespace Dmrg
#endif // KRONECKERDUMPER_H
