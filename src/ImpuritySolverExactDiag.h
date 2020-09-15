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
#include "InputNg.h"
#include "../../dmrgpp/src/Engine/InputCheck.h"
#include "LanczosSolver.h"
#include "MersenneTwister.h"
#include "BLAS.h"
#include "Matsubaras.h"

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
	typedef PsimagLite::InputNg<Dmrg::InputCheck> InputNgType;
	typedef PsimagLite::ParametersForSolver<RealType> SolverParametersType;
	typedef PsimagLite::LanczosSolver<SolverParametersType, SparseMatrixType, VectorComplexType>
	LanczosSolverType;
	typedef BasisType::LabeledOperatorType LabeledOperatorType;
	typedef PsimagLite::Matrix<ComplexOrRealType> MatrixType;
	typedef Matsubaras<RealType> MatsubarasType;

	ImpuritySolverExactDiag(const ParamsDmftSolverType& params,
	                        const ApplicationType& app)
	    : params_(params),
	      solverParams_(nullptr),
	      rng_(1234),
	      matsubaras_(params.ficticiousBeta, params.nMatsubaras),
	      gimp_(matsubaras_.total())
	{
		Dmrg::InputCheck inputCheck;
		InputNgType::Writeable ioW(params.gsTemplate, inputCheck);
		InputNgType::Readable io(ioW);
		io.read(hubbardU_, "hubbardU");
		io.readline(nup_, "TargetElectronsUp=");
		io.readline(ndown_, "TargetElectronsDown=");
		solverParams_ = new SolverParametersType(io, "Lanczos");
	}

	~ImpuritySolverExactDiag()
	{
		delete solverParams_;
		solverParams_ = nullptr;
	}

	// bathParams[0-nBath-1] ==> V ==> hoppings impurity --> bath
	// bathParams[nBath-...] ==> energies on each bath site
	void solve(const VectorRealType& bathParams)
	{
		ModelParamsType mp(bathParams);

		BasisType basis(mp.sites, nup_, ndown_);

		// setup model params ==> mp_

		SparseMatrixType matrix;
		setupHamiltonian(matrix, basis, mp);

		// Lanczos ==> gs AND Egs

		LanczosSolverType lanczos(matrix, *solverParams_);

		const SizeType n = matrix.rows();
		VectorComplexType initialVector(n);

		for (SizeType i = 0; i < n; ++i)
			initialVector[i] = rng_();

		RealType energy = 0;
		VectorComplexType gs(n);
		lanczos.computeOneState(energy, gs, initialVector, 1);

		// compute gimp
		computeGreenFunction(energy, gs, basis, mp);
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
				assert(hubbardU_.size() > i);
				s += hubbardU_[i] *
				        basis.isThereAnElectronAt(ket1,ket2,i,0,orb) * // SPIN_UP
				        basis.isThereAnElectronAt(ket1,ket2,i,1,orb); // SPIN_DOWN

				// Potential term
				RealType ne = (basis.getN(ket1,ket2,i,0,orb) +// SPIN_UP
				               basis.getN(ket1,ket2,i,1,orb));// SPIN_DOWN

				assert(mp.potentialV.size() > i);
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

	// <gs|c'(iwn-Hbar)^{-1}c|gs> + <gs|c(iwn+Hbar)^{-1}c'|gs>
	void computeGreenFunction(RealType energy,
	                          const VectorComplexType& gs,
	                          const BasisType& basis,
	                          const ModelParamsType& mp)
	{
		doType(0, energy, gs, basis, mp);

		doType(1, energy, gs, basis, mp);
	}

	// <gs|c' (iwn-Hbar)^{-1}c|gs> for what == 0; and the counterpart for what == 1
	void doType(SizeType what,
	            RealType energy,
	            const VectorComplexType& gs,
	            const BasisType& basis,
	            const ModelParamsType& mp)
	{
		SizeType center = mp.sites/2;
		SizeType spin = 0;

		MatrixType cAtCenter;
		SizeType nup = getElectrons(what, spin, nup_, 0);
		SizeType ndown = getElectrons(what, spin, ndown_, 1);

		BasisType basisDest(mp.sites, nup, ndown);
		LabeledOperatorType::Label label = (what == 0) ?
		            LabeledOperatorType::Label::OPERATOR_CDAGGER :
		            LabeledOperatorType::Label::OPERATOR_CDAGGER;
		setOperatorC(cAtCenter, basis, basisDest, label, center, spin);

		VectorComplexType opGs(basisDest.size());
		matrixVector(opGs, cAtCenter, gs, 'N');

		SparseMatrixType matrix;
		setupHamiltonian(matrix, basisDest, mp);
		MatrixType hmatrix = matrix.toDense();
		MatrixType oneOverMatrix = hmatrix;
		VectorComplexType correctionVector(basisDest.size());
		for (SizeType i = 0; i < matsubaras_.total(); ++i) {
			oneOverMatrix = hmatrix;
			const RealType wn = matsubaras_.omega(i);
			RealType sign = createOneOverMatrix(oneOverMatrix, matrix, wn, energy, what);
			matrixVector(correctionVector, oneOverMatrix, opGs, 'N');

			ComplexOrRealType sum = 0;
			for (SizeType site = 0; site < mp.sites; ++site) {
				MatrixType cAtSite;
				setOperatorC(cAtSite, basis, basisDest, label, site, spin);

				VectorComplexType opGs2(basisDest.size());
				matrixVector(opGs2, cAtSite, gs, 'N');
				ComplexOrRealType value = opGs2 * correctionVector;
				sum += value*sign;
			}

			if (what == 0)
				gimp_[i] = sum;
			else
				gimp_[i] += sum;
		}
	}

	//  -(-iwn + H - E0)  for what == 0
	//   +iwn + H - E0  for what == 1
	RealType createOneOverMatrix(MatrixType& oneOverMatrix,
	                             const SparseMatrixType& matrix,
	                             RealType energy,
	                             RealType wn,
	                             SizeType what) const
	{
		oneOverMatrix = matrix.toDense();
		RealType sign = (what == 0) ? -1 : 1;
		const SizeType n = oneOverMatrix.rows();
		for (SizeType i = 0; i < n; ++i) {
			oneOverMatrix(i, i) += ComplexType(-energy, sign*wn);
		}

		inverse(oneOverMatrix);
		// don't multiply oneOverMatrix by minus for what == 0, but keep track of sign instead
		return sign;
	}

	SizeType getElectrons(SizeType what,
	                      SizeType spin,
	                      SizeType electrons,
	                      SizeType upOrDown) const
	{
		if (spin != upOrDown) return electrons;
		return (what == 0) ? electrons - 1 : electrons + 1;
	}

	void setOperatorC(MatrixType& matrix,
	                  const BasisType& basisSrc,
	                  const BasisType& basisDest,
	                  LabeledOperatorType::Label label,
	                  SizeType site,
	                  SizeType spin) const
	{
		const SizeType hilbertSrc = basisSrc.size();
		const SizeType hilbertDest = basisDest.size();
		LabeledOperatorType lOperator(label);
		SizeType orb = 0;

		matrix.resize(hilbertDest, hilbertSrc);

		for (SizeType ispace = 0; ispace < hilbertSrc; ++ispace) {
			WordType ket1 = basisSrc(ispace, 0);
			WordType ket2 = basisSrc(ispace, 1);
			WordType bra = ket1;
			bool b = basisSrc.getBra(bra, ket1, ket2, lOperator, site, spin);
			if (!b) continue;
			SizeType index = basisDest.perfectIndex(bra, ket2);

			matrix(index, ispace) = basisDest.doSignGf(bra, ket2, site, spin, orb);
		}
	}

	void matrixVector(VectorComplexType& dest,
	                  const MatrixType& cAtCenter,
	                  const VectorComplexType& src,
	                  char trans) const
	{
		const ComplexType* gsPtr = &(src[0]);
		ComplexType* opGsPtr = &(dest[0]);

		//gemv (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
		psimag::BLAS::GEMV(trans,
		                   cAtCenter.rows(),
		                   cAtCenter.cols(),
		                   1,
		                   &(cAtCenter(0,0)),
		                   cAtCenter.rows(),
		                   gsPtr,
		                   1,
		                   0,
		                   opGsPtr,
		                   1);
	}

	const ParamsDmftSolverType& params_;
	SolverParametersType* solverParams_;
	PsimagLite::MersenneTwister rng_;
	MatsubarasType matsubaras_;
	VectorRealType hubbardU_;
	SizeType nup_;
	SizeType ndown_;
	VectorComplexType gimp_;
};
}
#endif // IMPURITYSOLVER_EXACTD_H
