#ifndef KRONECKERDUMPER_H
#define KRONECKERDUMPER_H
#include "Vector.h"
#include "TypeToString.h"
#include <fstream>
#include "../Version.h"
#include "Concurrency.h"
#include "ProgramGlobals.h"
#include "SymmetryElectronsSz.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class KroneckerDumper {

	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef SymmetryElectronsSz<RealType> SymmetryElectronsSzType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

public:

	struct ParamsForKroneckerDumper {
		ParamsForKroneckerDumper(bool enable = false,
		                         SizeType b = 0,
		                         SizeType e = 0,
		                         SizeType p = 6)
		    : enabled(enable), begin(b), end(e), precision(p)
		{}

		bool enabled;
		SizeType begin;
		SizeType end;
		SizeType precision;
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
		fout_.precision(p->precision);
		fout_<<"#KroneckerDumper for DMRG++ version "<<DMRGPP_VERSION<<"\n";
		fout_<<"#Instance="<<counter_<<"\n";
		fout_<<"#EncodingOfQuantumNumbers="<<(2*ProgramGlobals::maxElectronsOneSpin)<<"\n";

		printOneBasis("Left",lrs.left());
		printOneBasis("Right",lrs.right());

		fout_<<"#SuperBasisPermutation\n";
		fout_<<lrs.super().permutationVector();
		SizeType qtarget = lrs.super().qn(lrs.super().partition(m));
		PairSizeType etarget = getNupNdown(qtarget);
		fout_<<"#TargetElectronsUp="<<etarget.first<<"\n";
		fout_<<"#Target_ElectronsDown="<<etarget.second<<"\n";

		counter_++;
	}

	~KroneckerDumper()
	{
		fout_<<"#EOF\n";
		fout_.close();
	}

	void push(const SparseMatrixType& A, const SparseMatrixType& B, ComplexOrRealType val)
	{
		if (!enabled_) return;

		fout_<<"#START_AB_PAIR\n";
		fout_<<"link.value="<<val<<"\n";
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
		fout_<<matrix.row()<<" "<<matrix.col()<<"\n";
		for (SizeType i = 0; i < matrix.row(); ++i) {
			for (int k = matrix.getRowPtr(i); k < matrix.getRowPtr(i+1); ++k) {
				ComplexOrRealType value = matrix.getValue(k);
				if (PsimagLite::norm(value) == 0) continue;
				fout_<<i<<" "<<matrix.getCol(k)<<" "<<value<<"\n";
			}
		}
	}

	void printOneBasis(PsimagLite::String name, const BasisType& basis)
	{
		fout_<<"#" + name + "Basis\n";
		fout_<<"#Sites\n";
		fout_<<basis.block();
		fout_<<"#permutationVector\n";
		fout_<<basis.permutationVector();
		fout_<<"#ElectronsUp_ElectronsDown\n";
		fout_<<basis.size()<<"\n";
		for (SizeType i = 0; i < basis.size(); ++i) {
			SizeType q = basis.pseudoEffectiveNumber(i);
			PairSizeType nupDown = getNupNdown(q);
			fout_<<nupDown.first<<" "<<nupDown.second<<"\n";
		}

		fout_<<"#Electrons\n";
		fout_<<basis.electronsVector();
	}

	PairSizeType getNupNdown(SizeType q) const
	{
		VectorSizeType qns = SymmetryElectronsSzType::decodeQuantumNumber(q,2);
		SizeType electrons = qns[1];
		SizeType electronsUp = qns[0];
		assert(qns[1] >= qns[0]);
		SizeType electronsDown = electrons - qns[0];
		return PairSizeType(electronsUp,electronsDown);
	}

	static SizeType counter_;
	bool enabled_;
	std::ofstream fout_;
}; // class KroneckerDumpter

template<typename SparseMatrixType>
SizeType KroneckerDumper<SparseMatrixType>::counter_ = 0;

} // namespace Dmrg
#endif // KRONECKERDUMPER_H
