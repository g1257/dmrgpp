#ifndef OPERATORSTORAGE_H
#define OPERATORSTORAGE_H
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"
#include "Matrix.h"

// Selects storage for operators, they are blocked off diagonal,
// but blocks are now dense matrices, in the future we might make them
// MatrixDenseOrSparse type
// It also selects BlockDiagonalType for storage that we know is
// block diagonal, like the DMRG transformation matrix
namespace Dmrg {

template<typename ComplexOrRealType>
struct OperatorStorage {

	typedef Matrix<ComplexOrRealType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;
	typedef BlockDiagonalMatrixType BlockDiagonalType;
	typedef BlockOffDiagMatrixType BlockOffDiagType;
};
}
#endif // OPERATORSTORAGE_H
