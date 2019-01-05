#ifndef DMRG_APPLYHAMILTONIAN_H
#define DMRG_APPLYHAMILTONIAN_H

namespace Dmrg {

// ApplyHamiltonian, with Nsl = Not so local
template<typename ApplyOperatorLocalType,
         typename ModelType,
         typename LanczosSolverType>
class ApplyHamiltonian {

public:

	typedef typename ApplyOperatorLocalType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ApplyOperatorLocalType::OperatorType OperatorType;
	typedef typename ApplyOperatorLocalType::FermionSignType FermionSignType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename ModelType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename ApplyOperatorLocalType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::RealType RealType;

	ApplyHamiltonian(const LeftRightSuperType& lrs,
	                 const ModelType& model,
	                 const RealType& physicalTime)
	    : model_(model),
	      physicalTime_(physicalTime)
	{}

	//! FIXME: we need to make a fast version for when we're just
	//! figuring out where the (non-zero) partition is
	void operator()(VectorWithOffsetType& dest, const VectorWithOffsetType& phi) const
	{
		dest = phi;
		SizeType sectors = phi.sectors();
		for (SizeType ii = 0; ii < sectors; ++ii) {
			SizeType i = phi.sector(ii);
			VectorType r;
			internal_(r, phi, i);
			dest.setDataInSector(r, i);
		}
	}

private:

	void internal_(VectorType& r,
	               const VectorWithOffsetType& phi,
	               SizeType i0) const
	{
		SizeType p = applyOpLocal_.lrs().super().findPartitionNumber(phi.offset(i0));
		HamiltonianConnectionType hc(p,
		                             applyOpLocal_.lrs(),
		                             model_.geometry(),
		                             model_.modelLinks(),
		                             physicalTime_,
		                             0);
		MatrixForApplicationType lanczosHelper(model_, hc);

		SizeType total = phi.effectiveSize(i0);
		VectorType phi2(total);
		phi.extract(phi2, i0);
		r.resize(total);
		lanczosHelper.matrixVectorProduct(r, phi2);
	}

	ApplyHamiltonian(const ApplyHamiltonian&);

	ApplyHamiltonian& operator=(const ApplyHamiltonian&);

	ApplyOperatorLocalType applyOpLocal_;
	const ModelType& model_;
	const RealType& physicalTime_;
};
}
#endif // DMRG_APPLYHAMILTONIAN_H
