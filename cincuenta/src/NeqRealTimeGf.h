#ifndef NEQ_REAL_TIME_GF_H
#define NEQ_REAL_TIME_GF_H

#include <PsimagLite/Matrix.h>
#include <PsimagLite/Vector.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace Dmft {

/*! Real-time-only two-time Green's function on a uniform grid 0..N_t.
 *
 * This is the production GBEK representation.  Equilibrium Matsubara and
 * left-mixing data deliberately remain outside it: GBEK's propagation and
 * bath decomposition consume only G^R and G^<.
 */
template <typename ComplexOrRealType> class NeqRealTimeGf {

public:

	using RealType          = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using ComplexType       = std::complex<RealType>;
	using MatrixComplexType = PsimagLite::Matrix<ComplexType>;

	NeqRealTimeGf() = default;

	NeqRealTimeGf(SizeType nT, RealType dt)
	    : retarded(nT + 1, nT + 1)
	    , lesser(nT + 1, nT + 1)
	    , nT_(nT)
	    , dt_(dt)
	{ }

	SizeType nT() const { return nT_; }
	RealType dt() const { return dt_; }

	MatrixComplexType retarded; ///< G^R(t_n,t_j), causal for j <= n
	MatrixComplexType lesser; ///< G^<(t_n,t_j)

	void dump(const std::string& prefix) const
	{
		writeComponent(prefix + "-retarded", retarded, true);
		writeComponent(prefix + "-lesser", lesser, false);
	}

private:

	void writeComponent(const std::string&       filename,
	                    const MatrixComplexType& values,
	                    bool                     lowerTriangle) const
	{
		std::ofstream output(filename);
		if (!output)
			throw std::runtime_error("NeqRealTimeGf: cannot open " + filename);

		output << std::fixed << std::setprecision(10);
		for (SizeType n = 0; n <= nT_; ++n) {
			const SizeType last = lowerTriangle ? n : nT_;
			for (SizeType j = 0; j <= last; ++j) {
				const ComplexType value = values(n, j);
				output << n * dt_ << " " << j * dt_ << " " << value.real() << " "
				       << value.imag() << "\n";
			}
		}

		if (!output)
			throw std::runtime_error("NeqRealTimeGf: failed while writing " + filename);
	}

	SizeType nT_ = 0;
	RealType dt_ = 0;
};

} // namespace Dmft
#endif // NEQ_REAL_TIME_GF_H
