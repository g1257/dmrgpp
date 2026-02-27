#ifndef MATRIX_MARKET_HH
#define MATRIX_MARKET_HH
#include "CrsMatrix.h"
#include <string>

namespace Dmrg {

template <typename ComplexOrRealType> class MatrixMarket {

public:

	using VectorType       = std::vector<ComplexOrRealType>;
	using SparseMatrixType = PsimagLite::CrsMatrix<ComplexOrRealType>;

	MatrixMarket(const SparseMatrixType& sparse)
	    : sparse_(sparse)
	{ }

	void print(std::ostream& os) const
	{
		SizeType rows      = sparse_.rows();
		SizeType cols      = sparse_.cols();
		SizeType nonzeroes = sparse_.nonZeros();
		os << buildMMFormatHeader() << "\n";
		os << rows << " " << cols << " " << nonzeroes << "\n";
		for (SizeType i = 0; i < rows; ++i) {
			SizeType start = sparse_.getRowPtr(i);
			SizeType end   = sparse_.getRowPtr(i + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType                 col   = sparse_.getCol(k);
				const ComplexOrRealType& value = sparse_.getValue(k);
				// one-based output required by MM Format
				os << (i + 1) << " " << (col + 1) << " " << value << "\n";
			}
		}
	}

	std::string formatInfo()
	{
		return std::string(
		    "%============================================================================="
		    "====\n"
		    "%\n"
		    "% A sparse MxN matrix with L \n"
		    "% nonzeros in the Matrix Market format is represented as follows:\n"
		    "%\n"
		    "% +----------------------------------------------+\n"
		    "% |%%MatrixMarket matrix coordinate real general | <--- header line\n"
		    "% |%                                             | <--+\n"
		    "% |% comments                                    |    |-- 0 or more comment "
		    "lines\n"
		    "% |%                                             | <--+        \n"
		    "% |    M  N  L                                   | <--- rows, columns, "
		    "entries\n"
		    "% |    I1  J1  A(I1, J1)                         | <--+\n"
		    "% |    I2  J2  A(I2, J2)                         |    |\n"
		    "% |    I3  J3  A(I3, J3)                         |    |-- L lines\n"
		    "% |        . . .                                 |    |\n"
		    "% |    IL JL  A(IL, JL)                          | <--+\n"
		    "% +----------------------------------------------+   \n"
		    "%\n"
		    "% Indices are 1-based, i.e. A(1,1) is the first element.\n"
		    "%\n"
		    "%============================================================================="
		    "====\n");
	}

	static std::string buildMMFormatHeader()
	{
		return std::string("%%MatrixMarket matrix coordinate real general");
	}

	const SparseMatrixType& sparse_;
};
}
#endif // MATRIX_MARKET_HH
