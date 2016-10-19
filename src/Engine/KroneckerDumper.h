#ifndef KRONECKERDUMPER_H
#define KRONECKERDUMPER_H
#include "Vector.h"
#include "TypeToString.h"
#include <fstream>
#include "../Version.h"
#include "Concurrency.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class KroneckerDumper {

	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	struct ParamsForKroneckerDumper {
		ParamsForKroneckerDumper(bool enable = false,
		                         SizeType b = 0,
		                         SizeType e = 0)
		    : enabled(enable), begin(b), end(e)
		{}

		bool enabled;
		SizeType begin;
		SizeType end;
	}; // struct ParamsForKroneckerDumper

	KroneckerDumper(const ParamsForKroneckerDumper* p,
	                const LeftRightSuperType& lrs,
	                SizeType m)
	    : enabled_(p && p->enabled)
	{
		if (!enabled_) return;
		if (PsimagLite::Concurrency::npthreads > 1) {
			PsimagLite::String msg("KroneckerDumper cannot be run with Threads>1 ");
			throw PsimagLite::RuntimeError(msg + "because it's not thread-safe\n");
		}

		bool b = (p->end > 0 && counter_ >= p->end);
		if (counter_ < p->begin || b) {
			counter_++;
			enabled_ = false;
			return;
		}

		PsimagLite::String filename = "kroneckerDumper" + ttos(counter_) + ".txt";
		fout_.open(filename.c_str());
		fout_<<"#KroneckerDumper for DMRG++ version "<<DMRGPP_VERSION<<"\n";
		fout_<<"#Instance="<<counter_<<"\n";

		printOneBasis("Left",lrs.left());
		printOneBasis("Right",lrs.right());

		fout_<<"#TargetQuantumNumber="<<lrs.super().qn(lrs.super().partition(m))<<"\n";

		counter_++;
	}

	~KroneckerDumper()
	{
		fout_<<"#EOF\n";
		fout_.close();
	}

	void push(const SparseMatrixType& A, const SparseMatrixType& B)
	{
		if (!enabled_) return;

		fout_<<"#START_AB_PAIR\n";
		fout_<<"#A\n";
		printMatrix(A);
		fout_<<"#B\n";
		printMatrix(B);
		fout_<<"#END_AB_PAIR\n";
	}

	void push(bool option,const SparseMatrixType& hamiltonian)
	{
		if (!enabled_) return;
		if (option)
			fout_<<"#LeftHamiltonian\n";
		else
			fout_<<"#RightHamiltonian\n";
		printMatrix(hamiltonian);
	}

private:

	void printMatrix(const SparseMatrixType& matrix)
	{
		if (PRINTS_DENSE)
			fout_<<matrix.toDense();
		else
			fout_<<matrix;
	}

	void printOneBasis(PsimagLite::String name, const BasisType& basis)
	{
		fout_<<"#" + name + "Basis\n";
		fout_<<"#Sites\n";
		fout_<<basis.block();
		fout_<<"#permutationVector\n";
		fout_<<basis.permutationVector();
		fout_<<"#QuantumNumbers\n";
		for (SizeType i = 0; i < basis.size(); ++i)
			fout_<<basis.pseudoEffectiveNumber(i)<<"\n";
		fout_<<"#Electrons\n";
		fout_<<basis.electronsVector();
	}

	static const bool PRINTS_DENSE = true;
	static SizeType counter_;
	bool enabled_;
	bool printsDense_;
	std::ofstream fout_;
}; // class KroneckerDumpter

template<typename SparseMatrixType>
SizeType KroneckerDumper<SparseMatrixType>::counter_ = 0;

} // namespace Dmrg
#endif // KRONECKERDUMPER_H
