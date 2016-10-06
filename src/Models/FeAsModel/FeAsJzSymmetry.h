#ifndef FEASJZSYMMETRY_H
#define FEASJZSYMMETRY_H
#include "Vector.h"

namespace Dmrg {

template<typename HilbertBasisType, typename ComplexOrRealType>
class FeAsJzSymmetry {

public:

	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	FeAsJzSymmetry()
	{
		// write operator Jz in first basis
		// reorder Jz so that it is block diagonal --> permutation
		// diagonalize each block --> U', jzEigs
		// do U = U'*P^\dagger --> u_ and utranspose_
		// go through the blocks and determine the electrons
		// convertJzEigs(jzModifiedEigs_,jzEigs);
	}

	void setElectronsAndJz(VectorSizeType& electrons,
	                       VectorSizeType& electronsUp) const
	{
		if (!isEnabled_) return;
		electrons = electrons_;
		electronsUp = jzModifiedEigs_;
	}

	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType) const
	{
		if (basis.size() != electrons_.size())
			throw PsimagLite::RuntimeError("FeAsJzSymmetry: findElectrons");
		electrons = electrons_;
	}

	void jzReinterpret(MatrixType& cm) const
	{
		if (!isEnabled_) return;
		MatrixType tmp = utranspose_*cm;
		cm = tmp*u_;
	}

private:

	void convertJzEigs(VectorSizeType& electronselectronsUp,
	                   const VectorRealType& jzEigs) const
	{

	}

	bool isEnabled_;
	MatrixType u_;
	MatrixType utranspose_;
	VectorSizeType jzModifiedEigs_;
	VectorSizeType electrons_;
}; // class FeAsJzSymmetry
} // namespace Dmrg
#endif // FEASJZSYMMETRY_H
