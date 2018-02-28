#ifndef MINIAPPKRONECKER_H
#define MINIAPPKRONECKER_H
#include "Vector.h"
#include "IoSelector.h"
#include "VerySparseMatrix.h"
#include "Matrix.h"

namespace Dmrg {

template <typename ComplexOrRealType>
class MiniAppKronecker {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef VerySparseMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;

public:

	MiniAppKronecker(PsimagLite::String filename)
	{
		SparseMatrixType hamLeft(static_cast<SizeType>(0));
		SparseMatrixType hamRight(static_cast<SizeType>(0));
		PsimagLite::IoSelector::In io(filename);
		VectorSizeType pse;
		io.read(pse,"#SuperBasisPermutation");
		std::cerr<<"Read pse_.size="<<pse.size()<<"\n";
		pse_ = pse;
		for (SizeType i = 0; i < pse_.size(); ++i)
			pse_[pse[i]] = i;

		io.read(hamLeft, "#LeftHamiltonian");
		io.move(-20);
		std::cerr<<"Read H_L square, rank="<<hamLeft.rows()<<"\n";
		assert(isHermitian(hamLeft));

		io.read(hamRight, "#RightHamiltonian");
		io.move(-2);
		std::cerr<<"Read H_R square, rank="<<hamRight.rows()<<"\n";
		assert(isHermitian(hamRight));

		SizeType nLeft = hamLeft.rows();
		SizeType nRight = hamRight.cols();
		SizeType nSuper = nLeft*nRight;
		hamSuper_.resize(nSuper,nSuper);
		hamSuper_.setTo(0.0);
		buildHLeftAndRight(hamLeft,hamRight);
		assert(isHermitian(hamSuper_,true));


		SizeType counter = 0;
		while (!io.eof()) {
			SparseMatrixType Ahat(static_cast<SizeType>(0));
			SparseMatrixType B(static_cast<SizeType>(0));
			try {
				io.read(Ahat, "#Ahat" + ttos(counter));
			} catch (std::exception&) {
				break;
			}

			io.move(-20);
			io.read(B, "#B" + ttos(counter));
			io.move(-20);
			buildHconnection(Ahat,B);
			counter++;
		}

		assert(isHermitian(hamSuper_,true));
	}

	void printH(std::ostream& os) const
	{
		SizeType n = hamSuper_.n_row();
		os<<"P1\n";
		os<<n<<" "<<n<<"\n";
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				if (fabs(hamSuper_(i,j)) < 1e-6)
					os<<"0 ";
				else
					os<<"1 ";
			}

			os<<"\n";
		}
	}

private:

	void buildHconnection(const SparseMatrixType Ahat,
	                      const SparseMatrixType B)
	{
		SizeType nLeft = Ahat.rows();

		for (SizeType k = 0; k < Ahat.nonZeros(); ++k) {
			SizeType i = Ahat.getRow(k);
			SizeType ip = Ahat.getColumn(k);
			for (SizeType k2 = 0; k2 < B.nonZeros(); ++k2) {
				SizeType j = B.getRow(k2);
				SizeType jp = B.getColumn(k2);
				SizeType is = pack(i,j,nLeft);
				SizeType js = pack(ip,jp,nLeft);
				hamSuper_(is,js) += Ahat(i,ip)*B(j,jp);
			}
		}
	}

	void buildHLeftAndRight(const SparseMatrixType hamLeft,
	                        const SparseMatrixType hamRight)
	{
		SizeType nLeft = hamLeft.rows();
		SizeType nRight = hamRight.cols();

		// hamLeft
		for (SizeType k = 0; k < hamLeft.nonZeros(); ++k) {
			SizeType i = hamLeft.getRow(k);
			SizeType ip = hamLeft.getColumn(k);
			for (SizeType j = 0; j < nRight; ++j) {
				SizeType is = pack(i,j,nLeft);
				SizeType js = pack(ip,j,nLeft);
				hamSuper_(is,js) += hamLeft(i,ip);
			}
		}

		std::cerr<<"Done hamLeft\n";

		// hamRight
		for (SizeType k = 0; k < hamRight.nonZeros(); ++k) {
			SizeType i = hamRight.getRow(k);
			SizeType ip = hamRight.getColumn(k);
			for (SizeType j = 0; j < nLeft; ++j) {
				SizeType is = pack(j,i,nLeft);
				SizeType js = pack(j,ip,nLeft);
				hamSuper_(is,js) += hamRight(i,ip);
			}
		}

		std::cerr<<"Done hamRight\n";
	}

	SizeType pack(SizeType l, SizeType r, SizeType nLeft) const
	{
		assert(l + r*nLeft < pse_.size());
		return pse_[l + r*nLeft];
	}

	VectorSizeType pse_;
	MatrixType hamSuper_;
}; // class MiniAppKronecker
}

#endif // MINIAPPKRONECKER_H
