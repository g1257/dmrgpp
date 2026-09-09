#ifndef NEQ_HYBRIDIZATION_TARGET_H
#define NEQ_HYBRIDIZATION_TARGET_H

#include "NeqRealTimeGf.h"
#include "ParamsNeqDmftSolver.h"
#include <PsimagLite/PsimagLite.h>
#include <cmath>
#include <complex>

namespace Dmft {

/*!
 * \brief NeqHybridizationTarget
 * Bethe-lattice GBEK self-consistency target.
 *
 * The target supplied to the bath decomposition is the total hybridization
 * \f$\Lambda(t_n,t_j)=v(t_n)v(t_j)G_{imp}(t_n,t_j)\f$. It is not an
 * auxiliary propagator.
 */
template <typename ComplexOrRealType> class NeqHybridizationTarget {

public:

	using RealType       = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType    = std::complex<RealType>;
	using RealTimeGfType = NeqRealTimeGf<ComplexOrRealType>;
	using ParamsNeqType  = ParamsNeqDmftSolver<ComplexOrRealType>;

	explicit NeqHybridizationTarget(const ParamsNeqType& params)
	    : params_(params)
	    , vHop_(params.grid.v_hop)
	    , vHopSq_(vHop_ * vHop_)
	    , vHopFinal_(params.bandwidthFinal > RealType(0)
	                     ? RealType(0.25) * params.bandwidthFinal
	                     : vHop_)
	    , lambda_(params.nT, params.dt)
	{ }

	// Fill the complete causal row of the GBEK hybridization target,
	// Lambda(t_n,t_j) = v(t_n) v(t_j) Gimp(t_n,t_j).
	void updateTimeRow(int n, const RealTimeGfType& gimp)
	{
		const RealType vn = hoppingAt(n);
		for (int j = 0; j <= n; ++j) {
			const RealType vj      = hoppingAt(j);
			lambda_.retarded(n, j) = vn * vj * gimp.retarded(n, j);
			lambda_.lesser(n, j)   = vn * vj * gimp.lesser(n, j);
			lambda_.lesser(j, n)   = vj * vn * gimp.lesser(j, n);
		}
	}

	const RealTimeGfType& lambda() const { return lambda_; }

private:

	RealType hoppingAt(int n) const
	{
		const RealType tq = params_.quenchDuration;
		if (tq <= RealType(0) || params_.quenchShape == "step")
			return (n == 0) ? vHop_ : vHopFinal_;

		const RealType t = n * params_.dt;
		if (t >= tq)
			return vHopFinal_;

		RealType shape;
		if (params_.quenchShape == "tanh") {
			const RealType x = (t / tq - RealType(0.5)) * RealType(6);
			shape            = RealType(0.5) * (RealType(1) + std::tanh(x));
		} else {
			shape = RealType(0.5) * (RealType(1) - std::cos(M_PI * t / tq));
		}
		return vHop_ + (vHopFinal_ - vHop_) * shape;
	}

	const ParamsNeqType& params_;
	RealType             vHop_;
	RealType             vHopSq_;
	RealType             vHopFinal_;
	RealTimeGfType       lambda_;
};

} // namespace Dmft
#endif // NEQ_HYBRIDIZATION_TARGET_H
