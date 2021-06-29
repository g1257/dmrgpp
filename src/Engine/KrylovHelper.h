#ifndef KRYLOVHELPER_H
#define KRYLOVHELPER_H
#include "Vector.h"
#include "ProgressIndicator.h"

namespace Dmrg {

template<typename ActionType, typename TypeWrapperType>
class KrylovHelper {

public:

	typedef typename TypeWrapperType::MatrixComplexOrRealType MatrixComplexOrRealType;
	typedef typename TypeWrapperType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ActionType::VectorRealType VectorRealType;
	typedef typename TypeWrapperType::SolverParamsType SolverParamsType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;

	KrylovHelper(const SolverParamsType& params, SizeType firstRitz)
	    : params_(params), firstRitz_(firstRitz), progress_("KrylovHelper") {}

	template<typename SomeActionType>
	void calcR(VectorType& r,
	           const SomeActionType& whatRorI,
	           const MatrixComplexOrRealType& T,
	           const MatrixComplexOrRealType& V,
	           const VectorWithOffsetType& phi,
	           SizeType n2,
	           SizeType i0)
	{
		const bool krylovAbridge = !params_.options.isSet("KrylovNoAbridge");
		SizeType n3 = (krylovAbridge) ? 1 : n2;
		// ---------------------------------------------------
		// precompute values of calcVTimesPhi(kprime,v,phi,i0)
		// ---------------------------------------------------
		VectorType calcVTimesPhiArray(n3);
		for(SizeType kprime = 0; kprime < n3; ++kprime)
			calcVTimesPhiArray[kprime] = calcVTimesPhi(kprime, V, phi, i0);

		ComplexOrRealType sum2 = 0.0;
		for (SizeType k = firstRitz_; k < n2; ++k) {
			ComplexOrRealType sum = 0.0;
			for (SizeType kprime = 0; kprime < n3; ++kprime) {
				ComplexOrRealType tmp = PsimagLite::conj(T(kprime,k))*
				        calcVTimesPhiArray[kprime];
				sum += tmp;
				if (kprime > 0) sum2 += tmp;
			}

			r[k] = sum * whatRorI(k);
		}

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"Abridgment="<<sum2;
		if (krylovAbridge) msg<<" KrylovAbridge enabled";
		progress_.printline(msgg, std::cout);
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
	SizeType firstRitz_;
	PsimagLite::ProgressIndicator progress_;
};
}
#endif // KRYLOVHELPER_H
