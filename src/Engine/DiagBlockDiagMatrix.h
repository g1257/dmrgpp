#ifndef DIAGBLOCKDIAGMATRIX_H
#define DIAGBLOCKDIAGMATRIX_H
#include "EnforcePhase.h"

namespace Dmrg {

template<typename BlockDiagonalMatrixType>
class DiagBlockDiagMatrix {

	typedef typename BlockDiagonalMatrixType::BuildingBlockType BuildingBlockType;
	typedef typename BuildingBlockType::value_type ComplexOrRealType;
	typedef typename BlockDiagonalMatrixType::VectorRealType VectorRealType;

	class LoopForDiag {

		typedef PsimagLite::Concurrency ConcurrencyType;

	public:

		LoopForDiag(BlockDiagonalMatrixType& C1,
		            VectorRealType& eigs1,
		            char option1)
		    : C(C1),
		      eigs(eigs1),
		      option(option1),
		      eigsForGather(C.blocks()),
		      weights(C.blocks())
		{

			for (SizeType m=0;m<C.blocks();m++) {
				eigsForGather[m].resize(C.offsetsRows(m+1)-C.offsetsRows(m));
				weights[m] =  C.offsetsRows(m+1)-C.offsetsRows(m);
			}

			assert(C.rows() == C.cols());
			eigs.resize(C.rows());
		}

		SizeType tasks() const { return C.blocks(); }

		void doTask(SizeType taskNumber, SizeType)
		{
			assert(C.rows() == C.cols());
			SizeType m = taskNumber;
			VectorRealType eigsTmp;
			C.diagAndEnforcePhase(m, eigsTmp, option);
			for (SizeType j = C.offsetsRows(m); j < C.offsetsRows(m+1); ++j)
				eigsForGather[m][j-C.offsetsRows(m)] = eigsTmp[j-C.offsetsRows(m)];

		}

		void gather()
		{
			assert(C.rows() == C.cols());
			for (SizeType m = 0; m < C.blocks(); ++m) {
				for (SizeType j = C.offsetsRows(m);j < C.offsetsRows(m+1); ++j)
					eigs[j]=eigsForGather[m][j-C.offsetsRows(m)];
			}
		}

	private:

		BlockDiagonalMatrixType& C;
		VectorRealType& eigs;
		char option;
		typename PsimagLite::Vector<VectorRealType>::Type eigsForGather;
		typename PsimagLite::Vector<SizeType>::Type weights;
	};

public:

	// Parallel version of the diagonalization of a block diagonal matrix
	// Note: In reality, Parallelization is disabled here because a LAPACK call
	//        is needed and LAPACK is not necessarily thread safe.
	// This function is NOT called by useSvd
	static void diagonalise(BlockDiagonalMatrixType& C,
	                        VectorRealType& eigs,
	                        char option)
	{
		typedef PsimagLite::NoPthreadsNg<LoopForDiag> ParallelizerType;
		typedef PsimagLite::Concurrency ConcurrencyType;
		SizeType savedNpthreads = ConcurrencyType::npthreads;
		ConcurrencyType::npthreads = 1;
		ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
		                              PsimagLite::MPI::COMM_WORLD,
		                              false);

		LoopForDiag helper(C,eigs,option);

		threadObject.loopCreate(helper); // FIXME: needs weights

		helper.gather();

		ConcurrencyType::npthreads = savedNpthreads;
	}
}; // class DiagBlockDiagMatrix

} // namespace Dmrg
#endif // DIAGBLOCKDIAGMATRIX_H
