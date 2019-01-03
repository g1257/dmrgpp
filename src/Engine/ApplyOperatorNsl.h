#ifndef APPLYOPERATORNSL_H
#define APPLYOPERATORNSL_H

namespace Dmrg {

// ApplyOperatorNsl, with Nsl = Not so local
template<typename ApplyOperatorLocalType,
         typename ModelType,
         typename LanczosSolverType>
class ApplyOperatorNsl {

public:

	typedef typename ApplyOperatorLocalType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ApplyOperatorLocalType::OperatorType OperatorType;
	typedef typename ApplyOperatorLocalType::FermionSignType FermionSignType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename ModelType::HamiltonianConnectionType HamiltonianConnectionType;
	typedef typename ApplyOperatorLocalType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::RealType RealType;

	ApplyOperatorNsl(const LeftRightSuperType& lrs,
	                 bool withLegacyBug,
	                 const ModelType& model,
	                 const RealType& physicalTime)
	    : applyOpLocal_(lrs, withLegacyBug),
	      model_(model),
	      physicalTime_(physicalTime)
	{}

	//! FIXME: we need to make a fast version for when we're just
	//! figuring out where the (non-zero) partition is
	void operator()(VectorWithOffsetType& dest,
	                const VectorWithOffsetType& phi,
	                const OperatorType& AA,
	                const typename  ApplyOperatorLocalType::FermionSignType& fermionSign,
	                SizeType systemOrEnviron,
	                typename ApplyOperatorLocalType::BorderEnum corner) const
	{
		if (AA.category == OperatorType::CategoryEnum::REGULAR)
			return applyOpLocal_(dest, phi, AA, fermionSign, systemOrEnviron, corner);

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
		typename LanczosSolverType::MatrixType lanczosHelper(model_, hc);

		SizeType total = phi.effectiveSize(i0);
		VectorType phi2(total);
		phi.extract(phi2, i0);
		r.resize(total);
		lanczosHelper.matrixVectorProduct(r, phi2);
	}

	ApplyOperatorNsl(const ApplyOperatorNsl&);

	ApplyOperatorNsl& operator=(const ApplyOperatorNsl&);

	ApplyOperatorLocalType applyOpLocal_;
	const ModelType& model_;
	const RealType& physicalTime_;
};
}
#endif // APPLYOPERATORNSL_H
