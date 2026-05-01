/*
 */

#ifndef HUBBARDLANCZOS_H
#define HUBBARDLANCZOS_H

#include "BasisHubbardLanczos.h"
#include "BitManip.h"
#include "HubbardHelper.h"
#include "LanczosGlobals.h"
#include "LanczosModelBase.hpp"
#include "ParametersModelHubbard.h"
#include "TypeToString.h"

namespace LanczosPlusPlus {

template <typename ComplexOrRealType, typename GeometryType, typename InputType>
class HubbardOneOrbital : public LanczosModelBase<ComplexOrRealType, GeometryType, InputType> {

	enum
	{
		SPIN_UP   = LanczosGlobals::SPIN_UP,
		SPIN_DOWN = LanczosGlobals::SPIN_DOWN
	};

public:

	typedef typename PsimagLite::Real<ComplexOrRealType>::Type           RealType;
	typedef PsimagLite::Matrix<ComplexOrRealType>                        MatrixType;
	typedef LanczosModelBase<ComplexOrRealType, GeometryType, InputType> BaseType;
	typedef ParametersModelHubbard<RealType, InputType>                  ParametersModelType;
	typedef BasisHubbardLanczos<GeometryType>                            BasisType;
	typedef typename BasisType::PairIntType                              PairIntType;
	typedef typename BasisType::BaseType                                 BasisBaseType;
	typedef typename BasisType::WordType                                 WordType;
	typedef typename BaseType::VectorSizeType                            VectorSizeType;
	typedef typename BaseType::SparseMatrixType                          SparseMatrixType;
	typedef typename BaseType::VectorType                                VectorType;
	typedef typename BasisType::LabeledOperatorType                      LabeledOperatorType;
	typedef HubbardHelper<BaseType, BasisType, ParametersModelType>      HubbardHelperType;
	typedef typename HubbardHelperType::VectorRahulOperatorType VectorRahulOperatorType;

	static int const FERMION_SIGN = BasisType::BasisType::FERMION_SIGN;

	HubbardOneOrbital(SizeType nup, SizeType ndown, InputType& io, const GeometryType& geometry)
	    : mp_(io)
	    , geometry_(geometry)
	    , basis_(geometry, nup, ndown)
	    , helper_(geometry, mp_)
	{ }

	~HubbardOneOrbital() { BaseType::deleteGarbage(garbage_); }

	SizeType size() const { return basis_.size(); }

	SizeType orbitals(SizeType) const { return 1; }

	void setupHamiltonian(SparseMatrixType& matrix) const { setupHamiltonian(matrix, basis_); }

	//! Gf. related functions below:
	void setupHamiltonian(SparseMatrixType& matrix, const BasisBaseType& basis) const
	{
		helper_.setupHamiltonian(matrix, basis);
	}

	void matrixVectorProduct(VectorType& x, const VectorType& y) const
	{
		matrixVectorProduct(x, y, basis_);
	}

	void
	matrixVectorProduct(VectorType& x, VectorType const& y, const BasisBaseType& basis) const
	{
		helper_.matrixVectorProduct(x, y, basis);
	}

	bool hasNewParts(std::pair<SizeType, SizeType>&       newParts,
	                 const std::pair<SizeType, SizeType>& oldParts,
	                 const LabeledOperator&               lOperator,
	                 SizeType                             spin,
	                 SizeType) const
	{
		if (lOperator.id() == LabeledOperator::Label::OPERATOR_C
		    || lOperator.id() == LabeledOperator::Label::OPERATOR_CDAGGER)
			return hasNewPartsCorCdagger(newParts, oldParts, lOperator, spin);

		if (lOperator.id() == LabeledOperator::Label::OPERATOR_SPLUS
		    || lOperator.id() == LabeledOperator::Label::OPERATOR_SMINUS)
			return hasNewPartsSplusOrSminus(newParts, oldParts, lOperator, spin);

		if (lOperator.id() == LabeledOperator::Label::OPERATOR_SZ)
			return false;

		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += PsimagLite::String("hasNewParts: unsupported operator ");
		str += lOperator.toString() + "\n";
		throw std::runtime_error(str.c_str());
	}

	const GeometryType& geometry() const { return geometry_; }

	const BasisType& basis() const { return basis_; }

	PsimagLite::String name() const { return __FILE__; }

	BasisType* createBasis(SizeType nup, SizeType ndown) const
	{
		BasisType* ptr = new BasisType(geometry_, nup, ndown);
		garbage_.push_back(ptr);
		return ptr;
	}

	void print(std::ostream& os) const { os << mp_; }

	void printOperators(std::ostream& os) const
	{
		SizeType nup   = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		os << "#SectorSource 2 " << nup << " " << ndown << "\n";
		SizeType spin = SPIN_UP;
		for (SizeType site = 0; site < geometry_.numberOfSites(); ++site)
			printOperatorC(site, spin, os);
	}

private:

	void printOperatorC(SizeType site, SizeType spin, std::ostream& os) const
	{
		SizeType nup   = basis_.electrons(SPIN_UP);
		SizeType ndown = basis_.electrons(SPIN_DOWN);
		if (nup == 0) {
			os << "#Operator_c_" << spin << "_" << site << "\n";
			os << "#SectorDest 0\n"; // bogus
			os << "#Matrix\n";
			os << "0 0\n";
			return;
		}

		BasisType*     basis = createBasis(nup - 1, ndown);
		VectorSizeType opt(2, 0);
		opt[0] = site;
		opt[1] = spin;
		MatrixType matrix;
		setupOperator(matrix, *basis, "c", opt);
		os << "#Operator_c_" << spin << "_" << site << "\n";
		os << "#SectorDest 2 " << (nup - 1) << " " << ndown << "\n";
		os << "#Matrix\n";
		os << matrix;
	}

	void setupOperator(MatrixType&           matrix,
	                   const BasisBaseType&  basis,
	                   PsimagLite::String    operatorName,
	                   const VectorSizeType& operatorOptions) const
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

		if (operatorOptions.size() < 2) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "operator " + operatorName + " needs at least two options\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType site = operatorOptions[0];
		if (site >= nsite) {
			PsimagLite::String str(__FILE__);
			str += " " + ttos(__LINE__) + "\n";
			str += "site requested " + ttos(site);
			str += " but number of sites= " + ttos(nsite) + "\n";
			throw PsimagLite::RuntimeError(str);
		}

		SizeType spin = operatorOptions[1];
		matrix.resize(hilbertSrc, hilbertDest);
		matrix.setTo(0.0);
		SizeType orb = 0;

		for (SizeType ispace = 0; ispace < hilbertSrc; ispace++) {
			WordType ket1 = basis_(ispace, SPIN_UP);
			WordType ket2 = basis_(ispace, SPIN_DOWN);
			WordType bra  = ket1;
			// assumes OPERATOR_C
			bool b = basis.getBra(bra, ket1, ket2, lOperator, site, spin);
			if (!b)
				continue;
			SizeType index = basis.perfectIndex(bra, ket2);

			matrix(ispace, index) = basis.doSignGf(bra, ket2, site, spin, orb);
		}
	}

	bool hasNewPartsCorCdagger(std::pair<SizeType, SizeType>&       newParts,
	                           const std::pair<SizeType, SizeType>& oldParts,
	                           const LabeledOperatorType&           lOperator,
	                           SizeType                             spin) const
	{
		int newPart1 = oldParts.first;
		int newPart2 = oldParts.second;
		int c        = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_C) ? -1 : 1;
		if (spin == SPIN_UP)
			newPart1 += c;
		else
			newPart2 += c;

		if (newPart1 < 0 || newPart2 < 0)
			return false;
		SizeType nsite = geometry_.numberOfSites();
		if (SizeType(newPart1) > nsite || SizeType(newPart2) > nsite)
			return false;
		if (newPart1 == 0 && newPart2 == 0)
			return false;
		newParts.first  = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	bool hasNewPartsSplusOrSminus(std::pair<SizeType, SizeType>&       newParts,
	                              const std::pair<SizeType, SizeType>& oldParts,
	                              const LabeledOperatorType&           lOperator,
	                              SizeType) const
	{
		int newPart1 = oldParts.first;
		int newPart2 = oldParts.second;

		int c = (lOperator.id() == LabeledOperatorType::Label::OPERATOR_SPLUS) ? 1 : -1;
		newPart1 += c;
		newPart2 -= c;

		if (newPart1 < 0 || newPart2 < 0)
			return false;

		SizeType nsite = geometry_.numberOfSites();
		if (static_cast<SizeType>(newPart1) > nsite
		    || static_cast<SizeType>(newPart2) > nsite)
			return false;

		newParts.first  = SizeType(newPart1);
		newParts.second = SizeType(newPart2);
		return true;
	}

	const ParametersModelType                             mp_;
	const GeometryType&                                   geometry_;
	BasisType                                             basis_;
	mutable typename PsimagLite::Vector<BasisType*>::Type garbage_;
	HubbardHelperType                                     helper_;
}; // class HubbardOneOrbital
} // namespace LanczosPlusPlus
#endif
