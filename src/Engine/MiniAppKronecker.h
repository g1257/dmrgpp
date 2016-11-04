#ifndef MINIAPPKRONECKER_H
#define MINIAPPKRONECKER_H
#include "Vector.h"
#include "IoSimple.h"
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
	    : hamLeft_(static_cast<SizeType>(0)),
	      hamRight_(static_cast<SizeType>(0))
	{
		PsimagLite::IoSimple::In io(filename);
		io.advance("#LeftBasis");
		io.read(electrons_,"#Electrons");
		std::cerr<<"Read electrons_.size="<<electrons_.size()<<"\n";
		io.read(pse_,"#SuperBasisPermutation");
		std::cerr<<"Read pse_.size="<<pse_.size()<<"\n";
		io.readMatrix(hamLeft_,"#LeftHamiltonian");
		std::cerr<<"Read H_L square, rank="<<hamLeft_.rank()<<"\n";
		io.rewind();
		io.readMatrix(hamRight_,"#RightHamiltonian");
		std::cerr<<"Read H_R square, rank="<<hamRight_.rank()<<"\n";
		buildH();
	}

	void printH(std::ostream& os) const
	{
		SizeType n = hamSuper_.n_row();
		os<<"P1\n";
		os<<n<<" "<<n<<"\n";
		for (SizeType i = 0; i < n; ++i) {
			for (SizeType j = 0; j < n; ++j) {
				if (hamSuper_(i,j) == 0)
					os<<"0 ";
				else
					os<<"1 ";
			}

			os<<"\n";
		}
	}

private:

	void buildH()
	{
		SizeType nLeft = hamLeft_.rank();
		SizeType nRight = hamRight_.rank();
		SizeType nSuper = nLeft*nRight;
		hamSuper_.resize(nSuper,nSuper);
		hamSuper_.setTo(0.0);

		// hamLeft
		for (SizeType i = 0; i < nLeft; ++i) {
			for (SizeType j = 0; j < nRight; ++j) { // j = jprime
				SizeType is = i + j*nLeft;
				for (SizeType ip = 0; ip < nLeft; ++ip) { // iprime
					SizeType js = ip + j*nLeft;
					hamSuper_(is,js) = hamLeft_(i,ip);
				}
			}
		}
	}

	VectorSizeType electrons_;
	VectorSizeType pse_;
	SparseMatrixType hamLeft_;
	SparseMatrixType hamRight_;
	MatrixType hamSuper_;
}; // class MiniAppKronecker
}

#endif // MINIAPPKRONECKER_H
