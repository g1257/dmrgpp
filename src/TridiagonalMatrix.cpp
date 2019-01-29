#include "TridiagonalMatrix.h"

namespace PsimagLite {

template<>
void TridiagonalMatrix<double>::diag2(TridiagonalMatrix<double>::VectorRealType& eigs,
                                      SizeType nn) const
{
	char jobz = 'N';
	int n = nn;
	TridiagonalMatrix<double>::VectorType z(1, 0);
	int lz = 1;
	int lwork = 10;
	int liwork = 10;
	int info = 0;
	std::vector<int> iwork(liwork);
	TridiagonalMatrix<double>::VectorType work(lwork);
	TridiagonalMatrix<double>::VectorType e = b_;
	eigs = a_;

	psimag::LAPACK::dstedc_(&jobz,
	                        &n,
	                        &(eigs[0]),
	        &(e[1]),
	        &(z[0]),
	        &lz,
	        &(work[0]),
	        &lwork,
	        &(iwork[0]),
	        &liwork,
	        &info);

	if (info == 0) return;
	std::cerr<<"dstedc_ failed with info = "<<info<<"\n";
	throw RuntimeError("dstedc_ FAILED\n");
}

template<>
void TridiagonalMatrix<float>::diag2(TridiagonalMatrix<float>::VectorRealType& eigs,
                                     SizeType nn) const
{
	char jobz = 'N';
	int n = nn;
	TridiagonalMatrix<float>::VectorType z(1, 0);
	int lz = 1;
	int lwork = 10;
	int liwork = 10;
	int info = 0;
	std::vector<int> iwork(liwork);
	TridiagonalMatrix<float>::VectorType work(lwork);
	TridiagonalMatrix<float>::VectorType e = b_;
	eigs = a_;

	psimag::LAPACK::sstedc_(&jobz,
	                        &n,
	                        &(eigs[0]),
	        &(e[1]),
	        &(z[0]),
	        &lz,
	        &(work[0]),
	        &lwork,
	        &(iwork[0]),
	        &liwork,
	        &info);

	if (info == 0) return;
	std::cerr<<"sstedc_ failed with info = "<<info<<"\n";
	throw RuntimeError("dstedc_ FAILED\n");
}

}
