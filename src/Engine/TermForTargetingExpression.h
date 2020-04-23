#ifndef TERMFORTARGETINGEXPRESSION_H
#define TERMFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "AuxForTargetingExpression.h"
#include "OneOperatorSpec.h"

namespace Dmrg {

template<typename TargetingBaseType>
class TermForTargetingExpression {

public:

	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef AuxForTargetingExpression<TargetingBaseType> AuxiliaryType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename PsimagLite::Vector<OneOperatorSpecType*>::Type VectorOneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename ApplyOperatorExpressionType::BorderEnumType BorderEnumType;

	TermForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false), aux_(aux), factor_(1.0) {}

	TermForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	    : finalized_(false), aux_(aux), factor_(1.0), vStr_(1, str)
	{}

	TermForTargetingExpression& operator=(const TermForTargetingExpression& other)
	{
		finalized_ = other.finalized_;
		factor_ = other.factor_;
		vStr_ = other.vStr_;
		return *this;
	}

	void multiply(const TermForTargetingExpression& other)
	{
		const SizeType n = other.vStr_.size();
		for (SizeType i = 0; i < n; ++i)
			vStr_.push_back(other.vStr_[i]);
	}

	void finalize()
	{
		if (finalized_) return;

		SizeType n = vStr_.size();
		if (n == 0)
			err("AlgebraForTargetingExpression: Cannot finalize an empty object\n");

		const SizeType coo = getCurrentCoO();
		PsimagLite::String ket;
		SizeType sitesEqualToCoo = 0;
		VectorSizeType discardedTerms;

		VectorStringType newVstr;
		for (SizeType ii = 0; ii < n; ++ii) {
			const SizeType i = n -ii - 1; // read vector backwards
			PsimagLite::String tmp = vStr_[i];
			if (tmp[0] == '|') { // it's a vector
				if (ket != "")
					err("More than one ket found in " + toString() + "\n");
				ket = tmp;
				if (i + 1 != n)
					err("Vector is not at the end in " + toString() + "\n");

				newVstr.push_back(tmp);
				continue; // == first read
			}

			// it's a matrix or a scalar
			if (tmp[0] == '-' || tmp[0] == '+' || (tmp[0] >= 65 && tmp[0] <= 74)) {
				err("Scalars not supported yet\n");
			}

			// it's a matrix
			if (ii != n - 1) {
				newVstr.push_back(tmp);
				continue; // apply in order only IMPORTANT
				// the last is the ket
			}

			SiteSplitType siteSplit = OneOperatorSpecType::extractSiteIfAny(tmp);
			if (!siteSplit.hasSiteString)
				err("Each op must have a site\n");

			SizeType site = OneOperatorSpecType::strToNumberOfFail(siteSplit.siteString);
			if (site != coo) { // can only apply at center of orthogonality (coo)
				newVstr.push_back(tmp);
				continue;
			}

			if (sitesEqualToCoo >= 1)
				err("Sites cannot be repeated in term\n");
			++sitesEqualToCoo;
			discardedTerms.push_back(i);

			tmp = siteSplit.root;
			OneOperatorSpecType opspec(tmp);
			const PsimagLite::String destKet = tmp + "*" + ket;
			OperatorType* op = new OperatorType(aux_.model().naturalOperator(opspec.label,
			                                    0, // FIXME TODO SDHS Immm
			                                    opspec.dof));
			oneOperator(destKet, ket, *op, site);
			delete op;
			op = 0;
		}

		// discarded terms and inversion
		const SizeType nnew = newVstr.size();
		vStr_.clear();
		for (SizeType i = 0; i < nnew; ++i)
			vStr_.push_back(newVstr[nnew - i - 1]);

		finalized_ = true;
	}

	PsimagLite::String toString() const
	{
		PsimagLite::String s;
		const SizeType n = vStr_.size();
		if (n == 0)
			err("toString returns empty\n");

		for (SizeType i = 0; i < n - 1; ++i)
			s += vStr_[i] + "*";
		return s + vStr_[n - 1];
	}

private:

	void oneOperator(PsimagLite::String destKet,
	                 PsimagLite::String srcKet,
	                 const OperatorType& op,
	                 SizeType site)
	{
		SizeType currentCoO = getCurrentCoO();

		if (site == currentCoO) {
			const VectorWithOffsetType& srcVwo = aux_.getCurrentVectorConst(srcKet);
			PsimagLite::String internalName = aux_.createTemporaryVector(destKet);
			VectorWithOffsetType& destVwo = aux_.getCurrentVectorNonConst(internalName);
			applyInSitu(destVwo, srcVwo, site, op);
			return;
		}

		err("oneOperator should never reach here\n");
		// Fetch ket at coo site --> into vec
		// Apply op to vec ---> vec2
		// Bring vec2 back to current coo
	}

	// returns A|src1>
	void applyInSitu(VectorWithOffsetType& dest,
	                 const VectorWithOffsetType& src1,
	                 SizeType site,
	                 const OperatorType& A)
	{
		typename PsimagLite::Vector<bool>::Type oddElectrons;
		aux_.model().findOddElectronsOfOneSite(oddElectrons,site);
		FermionSign fs(aux_.lrs().left(), oddElectrons);
		bool b1 = (site == 1 &&
		           aux_.direction() == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);
		SizeType n = aux_.model().superGeometry().numberOfSites();
		assert(n > 2);
		bool b2 = (site == n - 2 &&
		           aux_.direction() == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		BorderEnumType border = (b1 || b2) ? BorderEnumType::BORDER_YES
		                                   : BorderEnumType::BORDER_NO;
		aux_.aoe().applyOpLocal()(dest, src1, A, fs, aux_.direction(), border);
	}

	SizeType getCurrentCoO() const
	{
		const LeftRightSuperType& lrs = aux_.lrs();
		const SizeType systemBlockSize = lrs.left().block().size();
		assert(systemBlockSize > 0);
		const int maxSystemSite = lrs.left().block()[systemBlockSize - 1];
		return (aux_.direction() == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? maxSystemSite
		                                                                        : maxSystemSite + 1;

	}

	 // ATTENTION: has assignment operator
	bool finalized_;
	const AuxiliaryType& aux_;
	ComplexOrRealType factor_;
	VectorStringType vStr_;
};
}
#endif // TERMFORTARGETINGEXPRESSION_H
