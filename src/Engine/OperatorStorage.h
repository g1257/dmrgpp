#ifndef OPERATORSTORAGE_H
#define OPERATORSTORAGE_H
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"
#include "Matrix.h"
#include "Io/IoNg.h"

// Selects storage for operators, they are blocked off diagonal,
// but blocks are now dense matrices, in the future we might make them
// MatrixDenseOrSparse type
// It also selects BlockDiagonalType for storage that we know is
// block diagonal, like the DMRG transformation matrix
namespace Dmrg {

template<typename ComplexOrRealType>
class OperatorStorage {

public:

	typedef ComplexOrRealType value_type;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;
	typedef BlockDiagonalMatrixType BlockDiagonalType;
	typedef BlockOffDiagMatrixType BlockOffDiagType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;

	OperatorStorage()
	{
		throw PsimagLite::RuntimeError("OperatorStorage::ctor(void)\n");
	}

	explicit OperatorStorage(const SparseMatrixType& src)
	{
		throw PsimagLite::RuntimeError("OperatorStorage::ctor(SparseMatrixType)\n");
	}

	void makeDiagonal(SizeType rows, ComplexOrRealType value = 1) // replace this by a ctor
	{
		throw PsimagLite::RuntimeError("OperatorStorage::makeDiagonal\n");
	}

	void read(PsimagLite::String label,
	          PsimagLite::IoNgSerializer& io)
	{
		throw PsimagLite::RuntimeError("OperatorStorage::read\n");
	}

	void write(PsimagLite::String label,
	          PsimagLite::IoNgSerializer& io,
	          PsimagLite::IoSerializer::WriteMode mode = PsimagLite::IoNgSerializer::NO_OVERWRITE)
	const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::write\n");
	}

	OperatorStorage operator+=(const OperatorStorage& other)
	{
		throw PsimagLite::RuntimeError("OperatorStorage::operator+=\n");
	}

	OperatorStorage operator*=(const ComplexOrRealType&)
	{
		throw PsimagLite::RuntimeError("OperatorStorage::operator*=\n");
	}

	void clear()
	{
		throw PsimagLite::RuntimeError("OperatorStorage::clear()\n");
	}

	bool checkValidity() const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::checkValidity\n");
	}

	void rotate(const PsimagLite::CrsMatrix<ComplexOrRealType>&,
	            const PsimagLite::CrsMatrix<ComplexOrRealType>&)
	{
		throw PsimagLite::RuntimeError("OperatorStorage::rotate\n");
	}

	MatrixType toDense() const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::toDense\n");
	}

	SparseMatrixType toCRS() const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::toCRS\n");
	}

	SizeType nonZeros() const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::nonZeros\n");
	}

	SizeType rows() const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::rows()\n");
	}

	SizeType cols() const
	{
		throw PsimagLite::RuntimeError("OperatorStorage::cols()\n");
	}
};

template<typename ComplexOrRealType>
OperatorStorage<ComplexOrRealType>
operator*(const typename OperatorStorage<ComplexOrRealType>::RealType& value,
          const OperatorStorage<ComplexOrRealType>& storage)
{
	throw PsimagLite::RuntimeError("OperatorStorage: operator*\n");
}

template<typename ComplexOrRealType>
OperatorStorage<ComplexOrRealType>
operator*(const OperatorStorage<ComplexOrRealType>&,
          const OperatorStorage<ComplexOrRealType>&)
{
	throw PsimagLite::RuntimeError("OperatorStorage: operator*\n");
}

template<typename T>
void fromCRS(OperatorStorage<T>& dest, const PsimagLite::CrsMatrix<T>& src)
{
	throw PsimagLite::RuntimeError("OperatorStorage: fromCRS\n");
}

template<typename T>
void toCRS(PsimagLite::CrsMatrix<T>& dest, const OperatorStorage<T>& src)
{
	dest = src.toCRS();
}

template<typename ComplexOrRealType>
void transposeConjugate(OperatorStorage<ComplexOrRealType>& dest,
                        const OperatorStorage<ComplexOrRealType>& src)
{
	err("OperatorStorage: transposeConjugate\n");
}

template<typename ComplexOrRealType>
void bcast(OperatorStorage<ComplexOrRealType>& dest)
{
	err("OperatorStorage: bcast\n");
}

template<typename ComplexOrRealType>
void crsMatrixToFullMatrix(PsimagLite::Matrix<ComplexOrRealType>& dest,
                           const OperatorStorage<ComplexOrRealType>& src)
{
	err("OperatorStorage: crsMatrixToFullMatrix\n");
}

template<typename ComplexOrRealType>
void fullMatrixToCrsMatrix(OperatorStorage<ComplexOrRealType>& dest,
                           const PsimagLite::Matrix<ComplexOrRealType>& src)
{
	err("OperatorStorage: fullMatrixToCrsMatrix\n");
}

template<typename ComplexOrRealType>
PsimagLite::Matrix<ComplexOrRealType> multiplyTc(const OperatorStorage<ComplexOrRealType>& src1,
                                                 const OperatorStorage<ComplexOrRealType>& src2)
{
	throw PsimagLite::RuntimeError("OperatorStorage: multiplyTc\n");
}

template<typename ComplexOrRealType>
bool isHermitian(const OperatorStorage<ComplexOrRealType>&)
{
	throw PsimagLite::RuntimeError("OperatorStorage: isHermitian\n");
}

template<typename ComplexOrRealType>
bool isAntiHermitian(const OperatorStorage<ComplexOrRealType>&)
{
	throw PsimagLite::RuntimeError("OperatorStorage: isHermitian\n");
}

template<typename ComplexOrRealType>
bool isTheIdentity(const OperatorStorage<ComplexOrRealType>&)
{
	throw PsimagLite::RuntimeError("OperatorStorage: isTheIdentity\n");
}

// See CrsMatrix.h line 734
template<typename T>
void externalProduct2(OperatorStorage<T>& B,
                      const OperatorStorage<T>& A,
                      SizeType nout,
                      const typename PsimagLite::Vector<typename PsimagLite::Real<T>::Type>::Type& signs,
                      bool order)
{
	throw PsimagLite::RuntimeError("OperatorStorage: externalProduct\n");
}
} // namespace PsimagLite

#endif // OPERATORSTORAGE_H
