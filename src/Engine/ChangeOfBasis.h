#ifndef DMRG_CHANGEOFBASIS_H
#define DMRG_CHANGEOFBASIS_H
namespace Dmrg {

template<typename SparseMatrixType, typename MatrixType>
class ChangeOfBasis {

public:

	void update(const SparseMatrixType& transform)
	{
		transform_ = transform;
		transposeConjugate(transformT_,transform);
	}

	void update(const MatrixType &v2)
	{
		SparseMatrixType v(v2);
		update(v);
	}

	void operator()(SparseMatrixType &v) const
	{
		SparseMatrixType tmp;
		multiply(tmp,v,transform_);
		multiply(v,transformT_,tmp);
	}

	static void changeBasis(SparseMatrixType &v,
	                        const SparseMatrixType& ftransform)
	{
		SparseMatrixType ftransformT;
		transposeConjugate(ftransformT,ftransform);
		SparseMatrixType tmp;
		multiply(tmp,v,ftransform);
		multiply(v,ftransformT,tmp);
	}

private:

	SparseMatrixType transform_;
	SparseMatrixType transformT_;
}; // class ChangeOfBasis
} // namespace Dmrg
#endif // DMRG_CHANGEOFBASIS_H
