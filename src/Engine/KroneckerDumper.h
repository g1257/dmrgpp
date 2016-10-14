#ifndef KRONECKERDUMPER_H
#define KRONECKERDUMPER_H
#include "Vector.h"

namespace Dmrg {

template<typename SparseMatrixType>
class KroneckerDumper {

	typedef typename PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	void push(const SparseMatrixType&,
	          const SparseMatrixType&,
	          const VectorBoolType&,
	          const VectorSizeType&,
	          SizeType start,
	          SizeType end) const
	{

	}
}; // class KroneckerDumpter

} // namespace Dmrg
#endif // KRONECKERDUMPER_H
