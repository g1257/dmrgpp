/*
 */

#ifndef LANCZOSPP_FERMION_SPINLESS_H
#define LANCZOSPP_FERMION_SPINLESS_H

#include "BasisFermionSpinless.hh"
#include "BitManip.h"
#include "LanczosGlobals.h"
#include "LanczosModelBase.hpp"
#include "Parallelizer2.h"
#include "ParametersFermionSpinless.hh"
#include "SparseRow.h"
#include "TypeToString.h"

namespace LanczosPlusPlus {

template <typename ComplexOrRealType, typename GeometryType, typename InputType>
class FermionSpinless : public LanczosModelBase<ComplexOrRealType, GeometryType, InputType> {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type           RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                        MatrixType;
	typedef LanczosModelBase<ComplexOrRealType, GeometryType, InputType> BaseType;
	typedef ParametersFermionSpinless<RealType, InputType>               ParametersModelType;
	typedef BasisFermionSpinless<GeometryType>                           BasisType;
	typedef typename BasisType::PairIntType                              PairIntType;
	typedef typename BasisType::BaseType                                 BasisBaseType;
	typedef typename BasisType::WordType                                 WordType;
	typedef typename BaseType::VectorSizeType                            VectorSizeType;
	typedef typename BaseType::SparseMatrixType                          SparseMatrixType;
	typedef typename BaseType::VectorType                                VectorType;
	typedef typename BasisType::LabeledOperatorType                      LabeledOperatorType;
	typedef PsimagLite::SparseRow<SparseMatrixType>                      SparseRowType;

	static int const FERMION_SIGN = BasisType::BasisType::FERMION_SIGN;

	FermionSpinless(SizeType ne, InputType& io, const GeometryType& geometry)
	    : mp_(io)
	    , geometry_(geometry)
	    , basis_(geometry, ne)
	    , hoppings_(geometry.numberOfSites(), geometry.numberOfSites())
	    , ninj_(geometry.numberOfSites(), geometry.numberOfSites())
	{
		const SizeType n = geometry_.numberOfSites();

		for (SizeType j = 0; j < n; ++j) {
			for (SizeType i = 0; i < n; ++i) {

				hoppings_(i, j) = geometry_(i, 0, j, 0, 0);
				ninj_(i, j)     = geometry_(i, 0, j, 0, 1);
			}
		}
	}

	~FermionSpinless() { BaseType::deleteGarbage(garbage_); }

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const { return 1; }

	void setupHamiltonian(SparseMatrixType& matrix) const { setupHamiltonian(matrix, basis_); }

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix, const BasisBaseType& basis) const
	{
		SizeType                                    hilbert = basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag, basis);

		SizeType nsite = geometry_.numberOfSites();

		matrix.resize(hilbert, hilbert);
		// Calculate off-diagonal elements AND store matrix
		SizeType nCounter = 0;
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			SparseRowType sparseRow;
			matrix.setRow(ispace, nCounter);
			WordType ket = basis(ispace, 0);
			// Save diagonal
			sparseRow.add(ispace, diag[ispace]);
			for (SizeType i = 0; i < nsite; ++i) {
				setHoppingTerm(sparseRow, ket, i, basis);
			}

			nCounter += sparseRow.finalize(matrix);
		}

		matrix.setRow(hilbert, nCounter);
	}

	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		matrixVectorProduct(x, y, basis_);
	}

	void
	matrixVectorProduct(VectorType& x, VectorType const& y, const BasisBaseType& basis) const
	{
		SizeType                                    hilbert = basis.size();
		typename PsimagLite::Vector<RealType>::Type diag(hilbert);
		calcDiagonalElements(diag, basis);
		for (SizeType ispace = 0; ispace < hilbert; ++ispace)
			x[ispace] += diag[ispace] * y[ispace];
		diag.clear();

		SizeType nsite = geometry_.numberOfSites();

		// Calculate off-diagonal elements AND store matrix
		auto lambda = [&basis, nsite, &x, &y, this](SizeType ispace, SizeType)
		{
			SparseRowType sparseRow;
			WordType      ket = basis(ispace, 0);
			for (SizeType i = 0; i < nsite; ++i) {
				setHoppingTerm(sparseRow, ket, i, basis);
			}

			x[ispace] += sparseRow.finalize(y);
		};

		PsimagLite::Parallelizer2<> parallelizer2(
		    PsimagLite::Concurrency::codeSectionParams);

		parallelizer2.parallelFor(0, hilbert, lambda);
	}

	bool hasNewParts(std::pair<SizeType, SizeType>&       newParts,
	                 const std::pair<SizeType, SizeType>& oldParts,
	                 const LabeledOperator&               lOperator,
	                 SizeType,
	                 SizeType) const
	{
		if (lOperator.id() == LabeledOperator::Label::OPERATOR_C
		    || lOperator.id() == LabeledOperator::Label::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts.first, oldParts.first, lOperator);

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += lOperator.toString() + "\n";
		throw std::runtime_error(str.c_str());
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType ne, SizeType) const
	{
		BasisType* ptr = new BasisType(geometry_, ne);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os << mp_; }

	void printOperators(std::ostream& os) const
	{
		SizeType ne = basis_.electrons();
		os << "#Electrons= " << ne << "\n";
		for (SizeType site = 0; site < geometry_.numberOfSites(); ++site)
			printOperatorC(site, os);
	}

private:

	void setHoppingTerm(SparseRowType&       sparseRow,
	                    const WordType&      ket,
	                    SizeType             i,
	                    const BasisBaseType& basis) const
	{
		WordType s1i = (ket & BasisType::bitmask(i));
		if (s1i > 0)
			s1i = 1;

		const SizeType nsite = geometry_.numberOfSites();

		// Hopping term
		for (SizeType j = 0; j < nsite; ++j) {
			const ComplexOrRealType& h = hoppings_(i, j);
			const bool hasHop = (PsimagLite::real(h) != 0 || PsimagLite::imag(h) != 0);
			WordType   s1j    = (ket & BasisType::bitmask(j));
			if (s1j > 0)
				s1j = 1;

			// Apply c^\dagger_j c_i
			if (hasHop && s1i == 1 && s1j == 0) {
				// apply i
				WordType bra1 = ket ^ BasisType::bitmask(i);
				RealType tmp2 = LanczosGlobals::doSign(ket, i)
				    * LanczosGlobals::doSign(bra1, j);

				// apply j
				bra1 = bra1 ^ BasisType::bitmask(j);

				SizeType          temp  = basis.perfectIndex(bra1, 0);
				ComplexOrRealType cTemp = h * tmp2; //*extraSign;
				assert(temp < basis.size());
				sparseRow.add(temp, cTemp);
			}
		}
	}

	void printOperatorC(SizeType site, std::ostream& os) const
	{
		SizeType ne = basis_.electrons();
		if (ne == 0) {
			os << "#Operator_c_" << site << "\n";
			os << "#SectorDest 0\n"; // bogus
			os << "#Matrix\n";
			os << "0 0\n";
			return;
		}

		BasisType* basis = createBasis(ne - 1, 0);
		MatrixType matrix;
		setupOperator(matrix, *basis, "c", site);
		os << "#Operator_c_" << site << "\n";
		os << "#SectorDest 2 " << (ne - 1) << "\n";
		os << "#Matrix\n";
		os << matrix;
	}

	void setupOperator(MatrixType&          matrix,
	                   const BasisBaseType& basis,
	                   PsimagLite::String   operatorName,
	                   SizeType             site) const
	{
		SizeType            hilbertDest = basis.size();
		SizeType            hilbertSrc  = basis_.size();
		SizeType            nsite       = geometry_.numberOfSites();
		LabeledOperatorType lOperator(LabeledOperatorType::Label::OPERATOR_C);
		if (operatorName != "c") {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " is unimplemented for this model\n";
			throw PsimagLite::RuntimeError(str);
		}

		if (site >= nsite) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "site requested " + ttos(site);
			str += " but number of sites= " + ttos(nsite) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		matrix.resize(hilbertSrc, hilbertDest);
		matrix.setTo(0.0);
		SizeType orb = 0;

		for (SizeType ispace = 0; ispace < hilbertSrc; ++ispace) {
			WordType ket = basis_(ispace, 0);
			WordType bra = ket;
			// assumes OPERATOR_C
			bool b = basis.getBra(bra, ket, 0, lOperator, site, 0);
			if (!b)
				continue;
			SizeType index = basis.perfectIndex(bra, 0);

			matrix(ispace, index) = basis.doSignGf(bra, 0, site, 0, orb);
		}
	}

	bool hasNewPartsCorCdagger(SizeType&                  newParts,
	                           SizeType                   oldParts,
	                           const LabeledOperatorType& lOperator) const
	{
		int newPart1 = oldParts;
		int c        = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C) ? -1 : 1;
		newPart1 += c;

		if (newPart1 < 0)
			return false;
		SizeType nsite = geometry_.numberOfSites();
		if (static_cast<SizeType>(newPart1) > nsite)
			return false;
		if (newPart1 == 0)
			return false;
		newParts = newPart1;
		return true;
	}

	void calcDiagonalElements(typename PsimagLite::Vector<RealType>::Type& diag,
	                          const BasisBaseType&                         basis) const
	{
		constexpr RealType zeroPointFive = 0.5;
		SizeType           hilbert       = basis.size();
		SizeType           nsite         = geometry_.numberOfSites();
		SizeType           orb           = 0;

		// Calculate diagonal elements
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			WordType          ket = basis(ispace, 0);
			ComplexOrRealType s   = 0;
			for (SizeType i = 0; i < nsite; ++i) {

				// ninj term
				RealType ne = basis.getN(ket, 0, i, 0, orb);

				for (SizeType j = 0; j < nsite; ++j) {
					ComplexOrRealType value = zeroPointFive * ninj_(i, j);
					if (PsimagLite::real(value) == 0
					    && PsimagLite::imag(value) == 0)
						continue;
					RealType tmp2 = basis.getN(ket, 0, j, 0, orb);
					s += value * ne * tmp2;
				}

				// Potential term
				RealType tmp = mp_.potentialV[i];
				s += tmp * ne;
			}

			assert(fabs(PsimagLite::imag(s)) < 1e-12);
			diag[ispace] = PsimagLite::real(s);
		}
	}

	const ParametersModelType                             mp_;
	const GeometryType&                                   geometry_;
	BasisType                                             basis_;
	MatrixType                                            hoppings_;
	MatrixType                                            ninj_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
}; // class FermionSpinless
} // namespace LanczosPlusPlus
#endif
