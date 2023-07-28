#ifndef OPERATORSTORAGE_H
#define OPERATORSTORAGE_H
#include "BlockDiagonalMatrix.h"
#include "BlockOffDiagMatrix.h"
#include "Io/IoNg.h"
#include "Matrix.h"

// Selects storage for operators,
// This can be just a CRS matrix or it can be the following.
// Blocked off diagonal matrix,
// but blocks are now dense matrices, in the future we might make them
// MatrixDenseOrSparse type
// It also selects BlockDiagonalType for storage that we know is
// block diagonal, like the DMRG transformation matrix
namespace Dmrg
{

template <typename ComplexOrRealType>
class OperatorStorage
{

public:

	typedef ComplexOrRealType value_type;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef BlockDiagonalMatrix<MatrixType> BlockDiagonalMatrixType;
	typedef BlockOffDiagMatrix<MatrixType> BlockOffDiagMatrixType;
	typedef BlockDiagonalMatrixType BlockDiagonalType;
	typedef BlockOffDiagMatrixType BlockOffDiagType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	enum class Type { CRS,
		NOT_READY };

	OperatorStorage(Type type = Type::CRS)
	    : type_(type)
	{
	}

	explicit OperatorStorage(const SparseMatrixType& src)
	    : type_(Type::CRS)
	    , crs_(src)
	{
	}

	void makeDiagonal(SizeType rows, ComplexOrRealType value = 1) // replace this by a ctor
	{
		if (!justCRS())
			throw PsimagLite::RuntimeError("OperatorStorage::makeDiagonal\n");

		crs_.makeDiagonal(rows, value);
	}

	void read(PsimagLite::String label,
	    PsimagLite::IoNgSerializer& io)
	{
		if (justCRS())
			return crs_.read(label, io);

		throw PsimagLite::RuntimeError("OperatorStorage::read\n");
	}

	void write(PsimagLite::String label,
	    PsimagLite::IoNgSerializer& io,
	    PsimagLite::IoSerializer::WriteMode mode = PsimagLite::IoNgSerializer::NO_OVERWRITE)
	    const
	{
		if (justCRS())
			return crs_.write(label, io, mode);

		throw PsimagLite::RuntimeError("OperatorStorage::write\n");
	}

	void overwrite(PsimagLite::String label,
	    PsimagLite::IoNgSerializer& io) const
	{
		if (justCRS())
			return crs_.overwrite(label, io);

		throw PsimagLite::RuntimeError("OperatorStorage::overwrite\n");
	}

	OperatorStorage operator+=(const OperatorStorage& other)
	{
		if (justCRS() && other.justCRS()) {
			crs_ += other.crs_;
			return *this;
		}

		throw PsimagLite::RuntimeError("OperatorStorage::operator+=\n");
	}

	OperatorStorage operator*=(const ComplexOrRealType& value)
	{
		if (justCRS()) {
			crs_ *= value;
			return *this;
		}

		throw PsimagLite::RuntimeError("OperatorStorage::operator*=\n");
	}

	void fromDense(const PsimagLite::Matrix<ComplexOrRealType>& m)
	{
		if (!justCRS())
			throw PsimagLite::RuntimeError("OperatorStorage::fromDense()\n");

		fullMatrixToCrsMatrix(crs_, m);
	}

	void clear()
	{
		if (justCRS())
			return crs_.clear();

		throw PsimagLite::RuntimeError("OperatorStorage::clear()\n");
	}

	void checkValidity() const
	{
		if (justCRS())
			return crs_.checkValidity();

		throw PsimagLite::RuntimeError("OperatorStorage::checkValidity\n");
	}

	void conjugate()
	{
		if (justCRS())
			return crs_.conjugate();

		throw PsimagLite::RuntimeError("OperatorStorage::conjugate\n");
	}

	void transpose()
	{
		if (!justCRS())
			throw PsimagLite::RuntimeError("OperatorStorage::transpose\n");

		// transpose conjugate
		SparseMatrixType copy = crs_;
		PsimagLite::transposeConjugate(crs_, copy);

		// conjugate again to end up transposing only
		crs_.conjugate();
	}

	void rotate(const PsimagLite::CrsMatrix<ComplexOrRealType>& left,
	    const PsimagLite::CrsMatrix<ComplexOrRealType>& right)
	{
		if (justCRS()) {
			SparseMatrixType tmp;
			multiply(tmp, crs_, right);
			multiply(crs_, left, tmp);
			return;
		}

		throw PsimagLite::RuntimeError("OperatorStorage::rotate\n");
	}

	MatrixType toDense() const
	{
		if (justCRS())
			return crs_.toDense();

		throw PsimagLite::RuntimeError("OperatorStorage::toDense\n");
	}

	const SparseMatrixType& getCRS() const
	{
		if (!justCRS())
			throw PsimagLite::RuntimeError("OperatorStorage::toCRS\n");
		return crs_;
	}

	// FIXME TODO DELETE THIS FUNCTION!!
	SparseMatrixType& getCRSNonConst()
	{
		if (!justCRS())
			throw PsimagLite::RuntimeError("OperatorStorage::toCRS\n");
		return crs_;
	}

	SizeType nonZeros() const
	{
		if (justCRS())
			return crs_.nonZeros();

		throw PsimagLite::RuntimeError("OperatorStorage::nonZeros\n");
	}

	SizeType rows() const
	{
		if (justCRS())
			return crs_.rows();

		throw PsimagLite::RuntimeError("OperatorStorage::rows()\n");
	}

	SizeType cols() const
	{
		if (justCRS())
			return crs_.cols();

		throw PsimagLite::RuntimeError("OperatorStorage::cols()\n");
	}

	bool justCRS() const { return (type_ == Type::CRS); }

	bool invalid() const { return (type_ == Type::NOT_READY); }

	friend void transposeConjugate(OperatorStorage& dest,
	    const OperatorStorage& src)
	{
		if (dest.justCRS() && src.justCRS())
			return transposeConjugate(dest.crs_, src.getCRS());

		err("OperatorStorage: transposeConjugate\n");
	}

	friend void fromCRS(OperatorStorage& dest,
	    const PsimagLite::CrsMatrix<ComplexOrRealType>& src)
	{
		if (dest.justCRS()) {
			dest.crs_ = src;
			return;
		}

		throw PsimagLite::RuntimeError("OperatorStorage: fromCRS\n");
	}

	friend void bcast(OperatorStorage& dest)
	{
		if (dest.justCRS())
			return bcast(dest.crs_);

		err("OperatorStorage: bcast\n");
	}

	// See CrsMatrix.h line 734
	friend void externalProduct2(OperatorStorage& B,
	    const OperatorStorage& A,
	    SizeType nout,
	    const VectorRealType& signs,
	    bool order,
	    const VectorSizeType& permutationFull)
	{
		if (B.justCRS() && A.justCRS())
			return externalProduct(B.crs_,
			    A.getCRS(),
			    nout,
			    signs,
			    order,
			    permutationFull);

		throw PsimagLite::RuntimeError("OperatorStorage: externalProduct\n");
	}

	friend void externalProduct2(OperatorStorage& C,
	    const OperatorStorage& A,
	    const OperatorStorage& B,
	    const VectorRealType& signs,
	    bool order,
	    const VectorSizeType& permutationFull)
	{
		if (A.justCRS() && B.justCRS() && C.justCRS())
			return externalProduct(C.crs_,
			    A.getCRS(),
			    B.getCRS(),
			    signs,
			    order,
			    permutationFull);

		throw PsimagLite::RuntimeError("OperatorStorage: externalProduct\n");
	}

	friend void fullMatrixToCrsMatrix(OperatorStorage& dest,
	    const PsimagLite::Matrix<ComplexOrRealType>& src)
	{
		if (dest.justCRS())
			return fullMatrixToCrsMatrix(dest.crs_, src);

		err("OperatorStorage: fullMatrixToCrsMatrix\n");
	}

private:

	Type type_;
	SparseMatrixType crs_;
};

template <typename ComplexOrRealType>
OperatorStorage<ComplexOrRealType>
operator*(const typename OperatorStorage<ComplexOrRealType>::RealType& value,
    const OperatorStorage<ComplexOrRealType>& storage)
{
	if (storage.justCRS())
		return storage.getCRS() * value;

	throw PsimagLite::RuntimeError("OperatorStorage: operator*\n");
}

template <typename ComplexOrRealType>
OperatorStorage<ComplexOrRealType>
operator*(const OperatorStorage<ComplexOrRealType>& a,
    const OperatorStorage<ComplexOrRealType>& b)
{
	if (a.justCRS() && b.justCRS())
		return OperatorStorage<ComplexOrRealType>(a.getCRS() * b.getCRS());

	throw PsimagLite::RuntimeError("OperatorStorage: operator*\n");
}

template <typename ComplexOrRealType>
void crsMatrixToFullMatrix(PsimagLite::Matrix<ComplexOrRealType>& dest,
    const OperatorStorage<ComplexOrRealType>& src)
{
	if (src.justCRS())
		return crsMatrixToFullMatrix(dest, src.getCRS());

	err("OperatorStorage: crsMatrixToFullMatrix\n");
}

template <typename ComplexOrRealType>
PsimagLite::Matrix<ComplexOrRealType> multiplyTc(const OperatorStorage<ComplexOrRealType>& src1,
    const OperatorStorage<ComplexOrRealType>& src2)
{
	if (src1.justCRS() && src2.justCRS())
		return multiplyTc(src1.getCRS(), src2.getCRS());

	throw PsimagLite::RuntimeError("OperatorStorage: multiplyTc\n");
}

template <typename ComplexOrRealType>
bool isHermitian(const OperatorStorage<ComplexOrRealType>& src)
{
	if (src.justCRS())
		return isHermitian(src.getCRS());

	throw PsimagLite::RuntimeError("OperatorStorage: isHermitian\n");
}

template <typename ComplexOrRealType>
bool isAntiHermitian(const OperatorStorage<ComplexOrRealType>& src)
{
	if (src.justCRS())
		return isAntiHermitian(src.getCRS());

	throw PsimagLite::RuntimeError("OperatorStorage: isHermitian\n");
}

template <typename ComplexOrRealType>
bool isTheIdentity(const OperatorStorage<ComplexOrRealType>& src)
{
	if (src.justCRS())
		return isTheIdentity(src.getCRS());

	throw PsimagLite::RuntimeError("OperatorStorage: isTheIdentity\n");
}

template <typename ComplexOrRealType>
bool isNonZeroMatrix(const OperatorStorage<ComplexOrRealType>& m)
{
	return (m.rows() > 0 && m.cols() > 0);
}

} // namespace PsimagLite

#endif // OPERATORSTORAGE_H
