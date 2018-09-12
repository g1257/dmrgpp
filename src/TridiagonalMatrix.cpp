#include "TridiagonalMatrix.h"

namespace PsimagLite {

template<>
void TridiagonalMatrix<double>::diag(TridiagonalMatrix<double>::VectorRealType& eigs,
                                     SizeType nn) const
{
	char jobz = 'N';
	int n = nn;
	TridiagonalMatrix<double>::VectorType z(1, 0);
	int lz = 1;
	int lwork = 10;
	int liwork = 10;
	int info = 0;
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
	        &liwork,
	        &liwork,
	        &info);

	if (info == 0) return;
	std::cerr<<"dstedc_ failed with info = "<<info<<"\n";
	throw RuntimeError("dstedc_ FAILED\n");
}

template<>
void TridiagonalMatrix<float>::diag(TridiagonalMatrix<float>::VectorRealType& eigs,
                                    SizeType nn) const
{
	char jobz = 'N';
	int n = nn;
	TridiagonalMatrix<float>::VectorType z(1, 0);
	int lz = 1;
	int lwork = 10;
	int liwork = 10;
	int info = 0;
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
	        &liwork,
	        &liwork,
	        &info);

	if (info == 0) return;
	std::cerr<<"sstedc_ failed with info = "<<info<<"\n";
	throw RuntimeError("dstedc_ FAILED\n");
}

/*template<>
void TridiagonalMatrix<std::complex<double> >::diag
(TridiagonalMatrix<std::complex<double> >::VectorRealType& eigs, SizeType nn) const
{
	char jobz = 'N';
	int n = nn;
	TridiagonalMatrix<std::complex<double> >::VectorType z(1, 0);
	int lz = 1;
	int lwork = 10;
	int liwork = 10;
	int info = 0;
	int lrwork = 10;
	TridiagonalMatrix<double>::VectorType rwork(lrwork);
	TridiagonalMatrix<double>::VectorType work(lwork);
	TridiagonalMatrix<double>::VectorType e = b_;
	eigs = a_;

	psimag::LAPACK::zstedc_(&jobz,
	                        &n,
	                        &(eigs[0]),
	        &(e[1]),
	        &(z[0]),
	        &lz,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &lrwork,
	        &liwork,
	        &liwork,
	        &info);

	if (info == 0) return;
	std::cerr<<"zstedc_ failed with info = "<<info<<"\n";
	throw RuntimeError("dstedc_ FAILED\n");
}

template<>
void TridiagonalMatrix<std::complex<float> >::diag
(TridiagonalMatrix<std::complex<float> >::VectorRealType& eigs, SizeType nn) const
{
	char jobz = 'N';
	int n = nn;
	TridiagonalMatrix<std::complex<float> >::VectorType z(1, 0);
	int lz = 1;
	int lwork = 10;
	int liwork = 10;
	int info = 0;
	int lrwork = 10;
	TridiagonalMatrix<float>::VectorType rwork(lrwork);
	TridiagonalMatrix<float>::VectorType work(lwork);
	TridiagonalMatrix<float>::VectorType e = b_;
	eigs = a_;

	psimag::LAPACK::cstedc_(&jobz,
	                        &n,
	                        &(eigs[0]),
	        &(e[1]),
	        &(z[0]),
	        &lz,
	        &(work[0]),
	        &lwork,
	        &(rwork[0]),
	        &lrwork,
	        &liwork,
	        &liwork,
	        &info);

	if (info == 0) return;
	std::cerr<<"cstedc_ failed with info = "<<info<<"\n";
	throw RuntimeError("dstedc_ FAILED\n");
}*/

}
