#ifndef BATCHEDGEMM_H
#define BATCHEDGEMM_H
#include "Vector.h"

namespace Dmrg {

template<typename InitKronType>
class BatchedGemm2 {

	typedef typename InitKronType::ArrayOfMatStructType ArrayOfMatStructType;
	typedef typename ArrayOfMatStructType::MatrixDenseOrSparseType MatrixDenseOrSparseType;
	typedef typename MatrixDenseOrSparseType::VectorType VectorType;
	typedef long int IntegerType;

public:

	BatchedGemm2(const InitKronType& initKron) : initKron_(initKron)
	{}

	bool enabled() const { return initKron_.batchedGemm(); }

	void matrixVector(VectorType&, const VectorType&) const
	{
		err("BatchedGemm: matrixVector not implemented yet, sorry\n");
	}

private:

	const InitKronType& initKron_;
	bool enabled_;
};
}
#endif // BATCHEDGEMM_H
