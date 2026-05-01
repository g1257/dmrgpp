#ifndef ONESECTOR_H
#define ONESECTOR_H
#include "BLAS.h"
#include "Matrix.h"
#include "Vector.h"

namespace LanczosPlusPlus {

template <typename RealType, typename InputType> class OneSector {

public:

	typedef PsimagLite::Vector<SizeType>::Type          VectorSizeType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Matrix<RealType>                MatrixType;

	OneSector(InputType& io)
	{
		io.read(sector_, "#SectorSource");
		io.read(eigs_, "#Eigenvalues");
		io.read(vecs_, "#Eigenvectors");
	}

	bool isSector(const VectorSizeType& jndVector) const { return (jndVector == sector_); }

	void info(std::ostream& os) const
	{
		os << "sector\n";
		os << sector_;
		os << "eigs.size()=" << eigs_.size() << "\n";
		os << "vecs=" << vecs_.n_row() << "x" << vecs_.n_col() << "\n";
	}

	SizeType size() const { return eigs_.size(); }

	void multiplyRight(MatrixType& x, const MatrixType& a) const
	{
		SizeType n = a.n_row();
		SizeType m = a.n_col();
		assert(vecs_.n_row() == m);
		assert(vecs_.n_col() == m);
		assert(x.n_row() == n);
		assert(x.n_col() == m);
		assert(m > 0 && n > 0);
		psimag::BLAS::GEMM(
		    'N', 'N', n, m, m, 1.0, &(a(0, 0)), n, &(vecs_(0, 0)), m, 0.0, &(x(0, 0)), n);
	}

	void multiplyLeft(MatrixType& x, const MatrixType& a) const
	{
		SizeType n = a.n_row();
		SizeType m = a.n_col();
		assert(vecs_.n_row() == n);
		psimag::BLAS::GEMM(
		    'C', 'N', n, m, n, 1.0, &(vecs_(0, 0)), n, &(a(0, 0)), n, 0.0, &(x(0, 0)), n);
	}

	const RealType& eig(SizeType i) const
	{
		assert(i < eigs_.size());
		return eigs_[i];
	}

private:

	VectorSizeType sector_;
	VectorRealType eigs_;
	MatrixType     vecs_;
};

} // namespace LanczosPlusPlus

#endif // ONESECTOR_H
