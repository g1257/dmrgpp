#ifndef FEASJZSYMMETRY_H
#define FEASJZSYMMETRY_H
#include "Vector.h"

namespace Dmrg {

template<typename HilbertBasisType, typename VectorOperatorType>
class FeAsJzSymmetry {

public:

	typedef typename VectorOperatorType::value_type OperatorType;
	typedef typename OperatorType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	FeAsJzSymmetry(bool isEnabled) :
	    isEnabled_(isEnabled), isSet_(false)
	{}

	void init(HilbertBasisType& natBasis,
	          VectorOperatorType& creationMatrix)
	{
		assert(!isSet_);
		// write operator Jz in first basis
		// reorder Jz so that it is block diagonal --> permutation
		// diagonalize each block --> U', jzEigs
		// do U = U'*P^\dagger --> u_ and utranspose_
		// go through the blocks and determine the electrons
		// convertJzEigs(jzModifiedEigs_,jzEigs);
		isSet_ = true;
	}

	void setElectronsAndJz(VectorSizeType& electrons,
	                       VectorSizeType& electronsUp) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		electrons = electrons_;
		electronsUp = jzModifiedEigs_;
	}

	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		if (basis.size() != electrons_.size())
			throw PsimagLite::RuntimeError("FeAsJzSymmetry: findElectrons");
		electrons = electrons_;
	}

	void jzReinterpret(MatrixType& cm) const
	{
		if (!isEnabled_) return;
		assert(isSet_);
		MatrixType tmp = utranspose_*cm;
		cm = tmp*u_;
	}

	bool isEnabled() const { return isEnabled_; }

	bool isSet() const { return isSet_; }

private:

	void convertJzEigs(VectorSizeType& electronselectronsUp,
	                   const VectorRealType& jzEigs) const
	{

	}

	bool isEnabled_;
	bool isSet_;
	MatrixType u_;
	MatrixType utranspose_;
	VectorSizeType jzModifiedEigs_;
	VectorSizeType electrons_;
}; // class FeAsJzSymmetry
} // namespace Dmrg
#endif // FEASJZSYMMETRY_H
