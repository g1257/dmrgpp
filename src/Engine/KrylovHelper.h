#ifndef KRYLOVHELPER_H
#define KRYLOVHELPER_H
#include "Vector.h"
#include "ProgressIndicator.h"

namespace Dmrg {

template<typename ActionType>
class KrylovHelper {

public:

	typedef typename ActionType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename ActionType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ActionType::VectorRealType VectorRealType;
	typedef typename ActionType::SolverParamsType SolverParamsType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;

	KrylovHelper(const SolverParamsType& params)
	    : params_(params), progress_("KrylovHelper") {}

	//	void calcR(VectorType& r,
	//	           const MatrixComplexOrRealType& T,
	//	           const MatrixComplexOrRealType& V,
	//	           const VectorWithOffsetType& phi,
	//	           RealType,
	//	           const VectorRealType& eigs,
	//	           SizeType timeIndex,
	//	           SizeType n2,
	//	           SizeType i0,
	//	           const TargetParamsType& tstStruct)
	//	{
	//		RealType timeDirection = tstStruct.timeDirection();

	//		for (SizeType k=0;k<n2;k++) {
	//			ComplexOrRealType sum = 0.0;
	//			for (SizeType kprime=0;kprime<n2;kprime++) {
	//				ComplexOrRealType tmpV = calcVTimesPhi(kprime,V,phi,i0);
	//				sum += PsimagLite::conj(T(kprime,k))*tmpV;
	//			}

	//			RealType tmp = (eigs[k]-E0_)*times_[timeIndex]*timeDirection;
	//			ComplexOrRealType c = 0.0;
	//			PsimagLite::expComplexOrReal(c,-tmp);
	//			r[k] = sum * c;
	//		}
	//	}

	void calcR(VectorType& r,
	           const ActionType& whatRorI,
	           const MatrixComplexOrRealType& T,
	           const MatrixComplexOrRealType& V,
	           const VectorWithOffsetType& phi,
	           SizeType n2,
	           SizeType i0)
	{
		bool krylovAbridge = (params_.options.find("KrylovNoAbridge") == PsimagLite::String::npos);
		SizeType n3 = (krylovAbridge) ? 1 : n2;
		// ---------------------------------------------------
		// precompute values of calcVTimesPhi(kprime,v,phi,i0)
		// ---------------------------------------------------
		VectorType calcVTimesPhiArray(n3);
		for(SizeType kprime = 0; kprime < n3; ++kprime)
			calcVTimesPhiArray[kprime] = calcVTimesPhi(kprime, V, phi, i0);

		ComplexOrRealType sum2 = 0.0;
		for (SizeType k = 0; k < n2; ++k) {
			ComplexOrRealType sum = 0.0;
			for (SizeType kprime = 0; kprime < n3; ++kprime) {
				ComplexOrRealType tmp = PsimagLite::conj(T(kprime,k))*
				        calcVTimesPhiArray[kprime];
				sum += tmp;
				if (kprime > 0) sum2 += tmp;
			}

			r[k] = sum * whatRorI(k);
		}

		PsimagLite::OstringStream msg;
		msg<<"Abridgment="<<sum2;
		if (krylovAbridge) msg<<" KrylovAbridge enabled";
		progress_.printline(msg, std::cout);
	}

	static ComplexOrRealType calcVTimesPhi(SizeType kprime,
	                                       const MatrixComplexOrRealType& V,
	                                       const VectorWithOffsetType& phi,
	                                       SizeType i0)
	{
		ComplexOrRealType ret = 0;
		SizeType total = phi.effectiveSize(i0);

		for (SizeType j = 0; j < total; ++j)
			ret += PsimagLite::conj(V(j,kprime))*phi.fastAccess(i0,j);
		return ret;
	}

private:

	const SolverParamsType& params_;
	PsimagLite::ProgressIndicator progress_;
};
}
#endif // KRYLOVHELPER_H
