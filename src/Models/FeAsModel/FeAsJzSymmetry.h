#ifndef FEASJZSYMMETRY_H
#define FEASJZSYMMETRY_H
#include "Vector.h"

namespace Dmrg {

template<typename HilbertBasisType, typename MatrixType>
class FeAsJzSymmetry {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	void setElectronsAndJzFor(VectorSizeType& electrons,
	                          VectorSizeType& electronsUp,
	                          SizeType ind) const
	{

	}

	void findElectrons(VectorSizeType& electrons,
	                   const HilbertBasisType& basis,
	                   SizeType site) const
	{}

	void jzReinterpret(MatrixType& cm) const
	{}
}; // class FeAsJzSymmetry
} // namespace Dmrg
#endif // FEASJZSYMMETRY_H
