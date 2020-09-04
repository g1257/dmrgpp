#ifndef IMPURITYSOLVER_EXACTD_H
#define IMPURITYSOLVER_EXACTD_H

#include "Vector.h"
#include "PsimagLite.h"
#include "ImpuritySolverBase.h"
#include "CrsMatrix.h"
#include "SparseRow.h"
#include "BitManip.h"
#include "ExactDiag/BasisExactDiag.h"
#include "ExactDiag/ModelParams.h"

namespace Dmft {

template<typename ParamsDmftSolverType>
class ImpuritySolverExactDiag : public ImpuritySolverBase<ParamsDmftSolverType> {

public:

	typedef typename ParamsDmftSolverType::ComplexOrRealType ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef std::complex<RealType> ComplexType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef typename PsimagLite::Vector<ComplexType>::Type VectorComplexType;
	typedef typename ImpuritySolverBase<ParamsDmftSolverType>::ApplicationType ApplicationType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType> SparseMatrixType;
	typedef PsimagLite::SparseRow<SparseMatrixType> SparseRowType;
	typedef long unsigned int WordType;
	typedef BasisExactDiag BasisType;
	typedef ModelParams<RealType> ModelParamsType;

	ImpuritySolverExactDiag(const ParamsDmftSolverType& params, const ApplicationType& app)
	    : params_(params)
	{}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{
		ModelParamsType mp(bathParams);
		SizeType nup = 0;
		SizeType ndown = 0;
		BasisType basis(mp.sites, nup, ndown);

		// setup model params ==> mp_

		SparseMatrixType matrix;
		setupHamiltonian(matrix, basis, mp);
		// matrix ==> dense
		PsimagLite::Matrix<ComplexOrRealType> mdense = matrix.toDense();

		VectorRealType eigs(mdense.rows());
		diag(mdense, eigs, 'V');

		// compute gimp
		err("Compute gimp unimplemented\n");
	}

	const VectorComplexType& gimp() const { return gimp_; }

private:

	void setupHamiltonian(SparseMatrixType& matrix,
	                      const BasisType& basis,
	                      const ModelParamsType& mp) const
	{
		SizeType hilbert = basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag, basis, mp);

		SizeType nsite = mp.sites;

		matrix.resize(hilbert, hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter=0;
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			SparseRowType sparseRow;
			matrix.setRow(ispace, nCounter);
			WordType ket1 = basis(ispace, 0); //SPIN_UP;
			WordType ket2 = basis(ispace, 1); //SPIN_DOWN;
			// Save diagonal
			sparseRow.add(ispace, diag[ispace]);
			for (SizeType i = 0; i < nsite; ++i) {
				setHoppingTerm(sparseRow, ket1, ket2, i, basis, mp);
			}

			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert,nCounter);
	}

	void calcDiagonalElements(VectorRealType& diag,
	                          const BasisType& basis,
	                          const ModelParamsType& mp) const
	{
		SizeType hilbert=basis.size();
		SizeType nsite = mp.sites;
		SizeType orb = 0;

		// Calculate diagonal elements
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			WordType ket1 = basis(ispace, 0); //SPIN_UP
			WordType ket2 = basis(ispace, 1); //SPIN_DOWN
			ComplexOrRealType s = 0;
			for (SizeType i = 0; i < nsite; ++i) {

				// Hubbard term U0
				s += mp.hubbardU[i] *
				        basis.isThereAnElectronAt(ket1,ket2,i,0,orb) * // SPIN_UP
				        basis.isThereAnElectronAt(ket1,ket2,i,1,orb); // SPIN_DOWN

				// Potential term
				RealType ne = (basis.getN(ket1,ket2,i,0,orb) +// SPIN_UP
				               basis.getN(ket1,ket2,i,1,orb));// SPIN_DOWN

				RealType tmp = mp.potentialV[i];
				if (tmp != 0) s += tmp * ne;
			}

			assert(fabs(PsimagLite::imag(s))<1e-12);
			diag[ispace] = PsimagLite::real(s);
		}
	}

	void setHoppingTerm(SparseRowType& sparseRow,
	                    const WordType& ket1,
	                    const WordType& ket2,
	                    SizeType i,
	                    const BasisType& basis,
	                    const ModelParamsType& mp) const
	{
		WordType s1i=(ket1 & BasisType::bitmask(i));
		if (s1i>0) s1i=1;
		WordType s2i=(ket2 & BasisType::bitmask(i));
		if (s2i>0) s2i=1;

		const SizeType nsite = mp.sites;

		// Hopping term
		for (SizeType j = 0; j < nsite; ++j) {
			const ComplexOrRealType& h = mp.hoppings(i, j);
			const bool hasHop = (PsimagLite::real(h) != 0 || PsimagLite::imag(h) != 0);
			WordType s1j= (ket1 & BasisType::bitmask(j));
			if (s1j>0) s1j=1;
			WordType s2j= (ket2 & BasisType::bitmask(j));
			if (s2j>0) s2j=1;

			// Apply c^\dagger_j c_i
			if (hasHop && s1i == 1 && s1j == 0) {
				// apply i
				WordType bra1 = ket1 ^ BasisType::bitmask(i);
				RealType tmp2 = doSign(ket1, i)*doSign(bra1, j);

				// apply j
				bra1 = bra1 ^ BasisType::bitmask(j);

				SizeType temp = basis.perfectIndex(bra1, ket2);
				//RealType extraSign = (s1j == 1) ? FERMION_SIGN : 1;
				ComplexOrRealType cTemp = h*tmp2; //*extraSign;
				//if (s1j == 1) cTemp = PsimagLite::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp, cTemp);
			}

			// Apply c^\dagger_j c_i DOWN
			if (hasHop && s2i == 1 && s2j == 0) {
				WordType bra2 = ket2 ^ BasisType::bitmask(i);
				RealType tmp2 = doSign(ket2, i)*doSign(bra2, j);

				bra2 = bra2 ^ BasisType::bitmask(j);

				SizeType temp = basis.perfectIndex(ket1, bra2);
				//RealType extraSign = (s2j == 1) ? FERMION_SIGN : 1;
				ComplexOrRealType cTemp = h*tmp2; //*extraSign;
				//if (s2j == 1) cTemp = PsimagLite::conj(cTemp);
				assert(temp<basis.size());
				sparseRow.add(temp, cTemp);
			}
		}
	}

	static int doSign(WordType a, SizeType i)
	{
		WordType mask = (1 << i) - 1;
		// Parity of single occupied between i and nsite-1
		return (PsimagLite::BitManip::count(a & mask) & 1) ? -1 : 1;
	}

	const ParamsDmftSolverType& params_;
	VectorComplexType gimp_;
};
}
#endif // IMPURITYSOLVER_EXACTD_H
