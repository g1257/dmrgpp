#ifndef PARAMS_MATSUBARA_GRID_H
#define PARAMS_MATSUBARA_GRID_H
#include "CincuentaInputCheck.h"
#include <PsimagLite/InputNg.h>
#include <PsimagLite/PsimagLite.h>
#include <PsimagLite/Vector.h>
#include <string>

namespace Dmft {

/*!
 * \brief ParamsMatsubaraGrid
 * The imaginary-time/Matsubara grid parameters that both the equilibrium
 * and non-equilibrium stages need: FicticiousBeta, ChemicalPotential,
 * Matsubaras, and (unless NeqAtomicLimit=1) LatticeGf.
 *
 * Split out of ParamsDmftSolver so the neq atomic-limit path -- which never
 * runs an equilibrium DMFT stage at all -- doesn't have to construct
 * equilibrium-bath-fit-only parameters (NumberOfBathPoints, ImpuritySolver,
 * FitOptions, MinParams*, ...) that would have no effect on that path.
 */
template <typename ComplexOrRealType> struct ParamsMatsubaraGrid {

	using RealType         = typename PsimagLite::Real<ComplexOrRealType>::Type;
	using InputNgType      = PsimagLite::InputNg<CincuentaInputCheck>;
	using VectorStringType = PsimagLite::Vector<PsimagLite::String>::Type;

	ParamsMatsubaraGrid(typename InputNgType::Readable& io, bool neqAtomicLimit)
	{
		io.readline(ficticiousBeta, "FicticiousBeta=");
		io.readline(mu, "ChemicalPotential=");
		io.readline(nMatsubaras, "Matsubaras=");

		// LatticeGf= sets the pre-quench (equilibrium) bandwidth; under
		// NeqAtomicLimit=1 there is no equilibrium bath at all, so it's not
		// just unused -- it's a stale key that would silently do nothing.
		// Reject it outright rather than requiring-and-ignoring it.
		if (neqAtomicLimit) {
			std::string tmp;
			bool        present = true;
			try {
				io.readline(tmp, "LatticeGf=");
			} catch (std::exception&) {
				present = false;
			}
			if (present)
				err("LatticeGf= has no effect when NeqAtomicLimit=1 (there is no "
				    "equilibrium bath to fit; the pre-quench hopping is exactly "
				    "zero by construction); remove it from the input\n");
		} else {
			io.readline(latticeGf, "LatticeGf=");
			v_hop = parseVhop(latticeGf);
		}
	}

	RealType    ficticiousBeta = 0;
	RealType    mu             = 0;
	SizeType    nMatsubaras    = 0;
	std::string latticeGf;

	/// v in Gramsch, Balzer, Eckstein, Kollar, PRB 88, 235106 (2013) Sec. VI /
	/// Eq. (14): the Bethe-lattice nearest-neighbor hopping (t_ij = v/sqrt(Z)
	/// in their notation), i.e. half the half-bandwidth (v = D/2 = W/4).
	/// Appears in the DMFT self-consistency Lambda(t,t') = v(t) G(t,t') v(t').
	/// Zero when NeqAtomicLimit=1 (no equilibrium bath at all).
	RealType v_hop = 0;

private:

	/*!
	 * \brief parseVhop
	 * Extract v = W/4 from "energy,semicircular,W" (see v_hop doc above).
	 * D = W/2 is the half-bandwidth; v = D/2 satisfies <epsilon^2> = v^2 for
	 * the semicircular DOS rho(eps) = sqrt(4v^2 - eps^2)/(2*pi*v^2) (GBEK
	 * Eq. 13). The neq Bethe-lattice self-consistency only supports this
	 * DOS shape, unlike the equilibrium LatticeGf machinery (LatticeGf.h),
	 * which also supports momentum-space dispersions and other densities of
	 * states -- so this check is strictly narrower than the equilibrium
	 * parse, not just a repeat of it.
	 *
	 * \param[in] latticeGf Lattice GF descriptor string, e.g. "energy,semicircular,W".
	 * \return v = W/4.
	 */
	static RealType parseVhop(const PsimagLite::String& latticeGf)
	{
		VectorStringType tokens;
		PsimagLite::split(tokens, latticeGf, ",");
		if (tokens.size() < 3 || tokens[0] != "energy" || tokens[1] != "semicircular")
			err("LatticeGf= must be 'energy,semicircular,W' for a non-equilibrium "
			    "run (the neq Bethe-lattice self-consistency only supports the "
			    "semicircular DOS); got: "
			    + latticeGf + "\n");
		return RealType(0.25) * PsimagLite::atof(tokens[2]);
	}
};

} // namespace Dmft
#endif // PARAMS_MATSUBARA_GRID_H
