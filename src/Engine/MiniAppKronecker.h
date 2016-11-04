#ifndef MINIAPPKRONECKER_H
#define MINIAPPKRONECKER_H
#include "Vector.h"
#include "IoSimple.h"
#include "VerySparseMatrix.h"

namespace Dmrg {

template <typename ComplexOrRealType>
class MiniAppKronecker {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef VerySparseMatrix<ComplexOrRealType> SparseMatrixType;

public:

	MiniAppKronecker(PsimagLite::String filename)
	    : hamLeft_(static_cast<SizeType>(0))
	{
		PsimagLite::IoSimple::In io(filename);
		io.advance("#LeftBasis");
		io.read(electrons_,"#Electrons");
		std::cerr<<"Read electrons_.size="<<electrons_.size()<<"\n";
		io.read(pse_,"#SuperBasisPermutation");
		std::cerr<<"Read pse_.size="<<pse_.size()<<"\n";
		io.readMatrix(hamLeft_,"#LeftHamiltonian");
		std::cerr<<"Read H_L square of rank="<<hamLeft_.rank();
		//buildH();
	}

	void check() const
	{

	}

private:

	void buildH()
	{
	}


	VectorSizeType electrons_;
	VectorSizeType pse_;
	SparseMatrixType hamLeft_;
}; // class MiniAppKronecker
}

#endif // MINIAPPKRONECKER_H
