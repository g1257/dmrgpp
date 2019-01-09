#ifndef ORACLECHEBYSHEV_H
#define ORACLECHEBYSHEV_H
#include "ScaledHamiltonian.h"

namespace Dmrg {

template<typename TargetingCommonType, typename TargetParamsType>
class OracleChebyshev {

public:

	typedef typename TargetingCommonType::ModelType ModelType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef typename TargetParamsType::RealType RealType;
	typedef typename TargetingCommonType::LanczosSolverType LanczosSolverType;
	typedef typename TargetingCommonType::ComplexOrRealType ComplexOrRealType;
	typedef typename TargetingCommonType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename LanczosSolverType::MatrixType MatrixLanczosType;
	typedef ScaledHamiltonian<MatrixLanczosType, TargetParamsType> ScaledHamiltonianType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename TargetingCommonType::FermionSignType FermionSignType;

	OracleChebyshev(const ModelType& model,
	                const LeftRightSuperType& lrs,
	                SizeType currentTime,
	                const TargetParamsType& tstStruct,
	                RealType E0)
	    : model_(model), lrs_(lrs), currentTime_(currentTime), tstStruct_(tstStruct), E0_(E0)
	{}

	void operator()(SizeType n,
	                const TargetingCommonType& common,
	                SizeType systemOrEnviron,
	                SizeType site,
	                OperatorType& A,
	                typename TargetingCommonType::BorderEnumType border)
	{
		VectorWithOffsetType p0;
		typename TargetingCommonType::ApplyOperatorType applyOpLocal(lrs_,
		                                                             common.withLegacyBugs());
		typename PsimagLite::Vector<bool>::Type signs;
		model_.findOddElectronsOfOneSite(signs, site);
		FermionSignType fs(lrs_.left(), signs);
		OperatorType Aprime = A;
		Aprime.dagger();
		applyOpLocal(p0,common.psi(),Aprime,fs,systemOrEnviron,border);

		VectorWithOffsetType p1;
		VectorType r;
		chebyIteration(r, p0, p0, true);
		p1.fromFull(r, lrs_.super());
		for (SizeType i = 0; i < n; ++i) {
			chebyIteration(r, p1, p0, false); // p2 = 2*H*p1 - p0;
			VectorWithOffsetType p2;
			p2.fromFull(r, lrs_.super());
			// <gs|c|p2>;
			ComplexOrRealType result = common.testRealWork(p2,
			                                               common.psi(),
			                                               systemOrEnviron,
			                                               site,
			                                               A,
			                                               border);
			std::cout<<"OracleChebyshev: <gs|H|p"<<(i+2)<<">= "<<result<<"\n";
			// prepare for next iteration
			p0 = p1;
			p1 = p2;
		}
	}

private:

	void chebyIteration(VectorType& r,
	                    const VectorWithOffsetType& p1,
	                    const VectorWithOffsetType& p0,
	                    bool firstOne) const
	{
		SizeType i0 = 0;
		for (SizeType ii = 0; ii < p1.sectors(); ++ii)
			i0 = p1.sector(ii);

		SizeType p = lrs_.super().findPartitionNumber(p1.offset(i0));
		typename ModelType::HamiltonianConnectionType hc(p,
		                                                 lrs_,
		                                                 model_.geometry(),
		                                                 ModelType::modelLinks(),
		                                                 currentTime_,
		                                                 0);
		MatrixLanczosType lanczosHelper(model_,
		                                hc);

		ScaledHamiltonianType lanczosHelper2(lanczosHelper, tstStruct_, E0_);

		SizeType total = p1.effectiveSize(i0);
		r.resize(total);
		VectorType x2(total);
		const RealType factor = (firstOne) ? 1.0 : 2.0;
		VectorWithOffsetType x = factor*p1;
		x.extract(x2, i0);
		lanczosHelper2.matrixVectorProduct(r, x2); // applying Hprime
		if (firstOne) return;

		VectorType phi2(total);
		p0.extract(phi2, i0);
 		r += (-1.0)*phi2;
	}

	const ModelType& model_;
	const LeftRightSuperType& lrs_;
	const SizeType& currentTime_;
	const TargetParamsType& tstStruct_;
	const RealType& E0_;
};
}
#endif // ORACLECHEBYSHEV_H
