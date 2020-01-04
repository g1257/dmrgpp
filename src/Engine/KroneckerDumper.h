#ifndef KRONECKERDUMPER_H
#define KRONECKERDUMPER_H
#include "Vector.h"
#include "TypeToString.h"
#include <fstream>
#include "../Version.h"
#include "Concurrency.h"
#include "ProgramGlobals.h"

namespace Dmrg {

template<typename LeftRightSuperType, typename SolverParamsType>
class KroneckerDumper {

	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename BasisType::QnType QnType;
	typedef std::pair<SizeType,SizeType> PairSizeType;

public:

	struct ParamsForKroneckerDumper {
		ParamsForKroneckerDumper(bool enable = false,
		                         SizeType b = 0,
		                         SizeType e = 0,
		                         SizeType p = 6,
		                         SizeType nOfQns_ = 0)
		    : enabled(enable), begin(b), end(e), precision(p),nOfQns(nOfQns_)
		{}

		bool enabled;
		SizeType begin;
		SizeType end;
		SizeType precision;
		SizeType nOfQns;
	}; // struct ParamsForKroneckerDumper

	KroneckerDumper(const SolverParamsType& params,
	                const LeftRightSuperType& lrs,
	                ProgramGlobals::DirectionEnum dir)
	    : enabled_(false),pairCount_(0),disable_(false)
	{
		if (dir == ProgramGlobals::DirectionEnum::INFINITE)
			return;

		enabled_ = params.options.isSet("KroneckerDumper");
		if (!enabled_) return;
		ParamsForKroneckerDumper p(enabled_,
		                           params.dumperBegin,
		                           params.dumperEnd,
		                           params.precision);

		bool b = (p.end > 0 && counter_ >= p.end);
		if (counter_ < p.begin || b) {
			++counter_;
			enabled_ = false;
			return;
		}

		if (p.nOfQns == 0) {
			PsimagLite::String msg("KroneckerDumper::ctor(): internal error ");
			throw PsimagLite::RuntimeError(msg + "nOfQns\n");
		}

		ConcurrencyType::mutexInit(&mutex_);

		PsimagLite::String filename = "kroneckerDumper" + ttos(counter_) + ".txt";
		fout_.open(filename.c_str());
		fout_.precision(p.precision);
		fout_<<"KroneckerDumper for DMRG++ version "<<DMRGPP_VERSION<<"\n";
		fout_<<"Instance="<<counter_<<"\n";
		fout_<<"EncodingOfQuantumNumbers="<<(2*ProgramGlobals::maxElectronsOneSpin)<<"\n";

		printOneBasis("Left",lrs.left(),p.nOfQns);
		printOneBasis("Right",lrs.right(),p.nOfQns);

		fout_<<"SuperBasisPermutation\n";
		fout_<<lrs.super().permutationVector();
		//QnType qtarget = lrs.super().qnEx(m);
		//fout_<<qtarget<<"\n";

		signs_ = lrs.left().signs();
		++counter_;
	}

	~KroneckerDumper()
	{
		if (!enabled_) return;

		fout_<<"EOF\n";
		fout_.close();

		ConcurrencyType::mutexDestroy(&mutex_);
	}

	void push(const SparseMatrixType& A,
	          const SparseMatrixType& B,
	          ComplexOrRealType val,
	          ProgramGlobals::FermionOrBosonEnum bosonOrFermion,
	          const VectorType& y)
	{
		if (!enabled_) return;
		if (disable_) return;

		assert(y_.size() > 0);
		if (notFirstVector(y)) {
			disable_ = true;
			return;
		}

		ConcurrencyType::mutexLock(&mutex_);
		fout_<<"START_AB_PAIR\n";
		fout_<<"link.value="<<val<<"\n";
		fout_<<"A"<<pairCount_<<"\n";
		printMatrix(A);
		fout_<<"Ahat"<<pairCount_<<"\n";
		SparseMatrixType Ahat;
		calculateAhat(Ahat,A,val,bosonOrFermion);
		printMatrix(Ahat);
		fout_<<"B"<<pairCount_<<"\n";
		printMatrix(B);
		fout_<<"END_AB_PAIR\n";
		pairCount_++;
		ConcurrencyType::mutexUnlock(&mutex_);
	}

	void push(bool option,
	          const SparseMatrixType& hamiltonian,
	          const VectorType& y)
	{
		if (!enabled_) return;
		if (disable_) return;

		if (y_.size() == 0) y_ = y;
		if (notFirstVector(y)) {
			disable_ = true;
			return;
		}

		if (option)
			fout_<<"LeftHamiltonian\n";
		else
			fout_<<"RightHamiltonian\n";
		printMatrix(hamiltonian);
	}

private:

	void printMatrix(const SparseMatrixType& matrix)
	{
		fout_<<matrix.rows()<<" "<<matrix.cols()<<"\n";
		for (SizeType i = 0; i < matrix.rows(); ++i) {
			for (int k = matrix.getRowPtr(i); k < matrix.getRowPtr(i+1); ++k) {
				ComplexOrRealType value = matrix.getValue(k);
				if (PsimagLite::norm(value) == 0) continue;
				fout_<<i<<" "<<matrix.getCol(k)<<" "<<value<<"\n";
			}
		}
	}

	void printOneBasis(PsimagLite::String name,
	                   const BasisType& basis,
	                   SizeType nOfQns)
	{
		fout_<<"" + name + "Basis\n";
		fout_<<"Sites\n";
		fout_<<basis.block();
		fout_<<"permutationVector\n";
		fout_<<basis.permutationVector();
		fout_<<"ElectronsUp_ElectronsDown\n";
		fout_<<basis.size()<<"\n";
		SizeType npart = basis.partition() - 1;
		for (SizeType i = 0; i < npart; ++i)
			fout_<<basis.pseudoQnToString(i)<<"\n";

		fout_<<"Signs\n";
		fout_<<basis.signs();
	}

	PairSizeType getNupNdown(QnType q) const
	{
		assert(q.other.size() >= 1);
		return PairSizeType(q.other[0], q.electrons - q.other[0]);
	}

	// Ahat(ia,ja) = (-1)^e_L(ia) A(ia,ja)*value
	void calculateAhat(SparseMatrixType& Ahat,
	                   const SparseMatrixType& A,
	                   ComplexOrRealType val,
	                   ProgramGlobals::FermionOrBosonEnum bosonOrFermion) const
	{
		Ahat = A;
		SizeType rows = Ahat.rows();
		SizeType counter = 0;
		for (SizeType i = 0; i < rows; ++i) {
			RealType sign = (bosonOrFermion == ProgramGlobals::FermionOrBosonEnum::FERMION &&
			                 signs_[i]) ? -1.0 : 1.0;
			for (int k = Ahat.getRowPtr(i); k < Ahat.getRowPtr(i+1); ++k) {
				ComplexOrRealType tmp = Ahat.getValue(k)*sign*val;
				Ahat.setValues(counter++, tmp);
			}
		}
	}

	bool notFirstVector(const VectorType& y) const
	{
		VectorType ydiff;
		ydiff <= y - y_;
		return (PsimagLite::norm(ydiff) > 1e-8);
	}

	static SizeType counter_;
	bool enabled_;
	SizeType pairCount_;
	bool disable_;
	VectorType y_;
	std::ofstream fout_;
	VectorBoolType signs_;
	ConcurrencyType::MutexType mutex_;
}; // class KroneckerDumpter

template<typename T1, typename T2>
SizeType KroneckerDumper<T1, T2>::counter_ = 0;

} // namespace Dmrg
#endif // KRONECKERDUMPER_H
