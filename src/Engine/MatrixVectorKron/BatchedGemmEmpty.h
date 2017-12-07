#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;

public:

	BatchedGemm(const InitKronType&) {}

	bool enabled() const { return false; }

	void matrixVector(VectorType&, const VectorType&) const
	{}
};
}
#endif // BATCHEDGEMM_H
