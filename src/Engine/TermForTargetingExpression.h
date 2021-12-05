#ifndef TERMFORTARGETINGEXPRESSION_H
#define TERMFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "AuxForTargetingExpression.h"
#include "OneOperatorSpec.h"
#include "NonLocalForTargetingExpression.h"

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
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef AuxForTargetingExpression<TargetingBaseType> AuxiliaryType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename PsimagLite::Vector<OneOperatorSpecType*>::Type VectorOneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename ApplyOperatorExpressionType::BorderEnumType BorderEnumType;
	typedef NonLocalForTargetingExpression<TargetingBaseType> NonLocalForTargetingExpressionType;

	TermForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false), aux_(aux), factor_(1.0), nonLocal_(aux) {}

	TermForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	    : finalized_(false), aux_(aux), factor_(1.0), vStr_(1, str), nonLocal_(aux)
	{}

	TermForTargetingExpression& operator=(const TermForTargetingExpression& other) = delete;

	TermForTargetingExpression(const TermForTargetingExpression& other) = delete;

	void assignAndDestroy(TermForTargetingExpression& other)
	{
		finalized_ = other.finalized_;
		factor_ = other.factor_;
		strFactor_ = other.strFactor_;
		vStr_ = other.vStr_;
	}

	void multiply(const TermForTargetingExpression& other)
	{
		const SizeType n = other.vStr_.size();
		for (SizeType i = 0; i < n; ++i)
			vStr_.push_back(other.vStr_[i]);
	}

	void multiply(ComplexOrRealType val)
	{
		factor_ *= val;

		return (PsimagLite::IsComplexNumber<ComplexOrRealType>::True) ? finalMultImag()
		                                                              : finalMultReal();
	}

	void finalize()
	{
		if (finalized_) return;

		SizeType n = vStr_.size();
		if (n == 0)
			err("AlgebraForTargetingExpression: Cannot finalize an empty object\n");

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
			if (tmp[0] == '-' || tmp[0] == '+' || tmp[0] == '.' || (tmp[0] >= 48 && tmp[0] <= 57)) {
				err("Scalars not supported yet\n");
			}

			// it's a matrix
			if (ii != 1) {
				newVstr.push_back(tmp);
				continue; // apply in order only IMPORTANT
				// the last is the ket
			}

			SiteSplitType siteSplit = OneOperatorSpecType::extractSiteIfAny(tmp);
			if (isGlobalOperator(tmp)) {
				nonLocal_.timeEvolve(siteSplit, ket, tmp, getCurrentCoO());
				continue;
			}

			if (!siteSplit.hasSiteString)
				err("This op must have a site\n");

			SizeType site = OneOperatorSpecType::strToNumberOrFail(siteSplit.siteString);
			// can only apply at center of orthogonality (coo)
			if (!siteCanBeApplied(site)) {
				newVstr.push_back(tmp);
				continue;
			}

			// for now we delay application of more than one site
			// even if possible at the current coo
			if (sitesEqualToCoo > 0) {
				newVstr.push_back(tmp);
				continue;
			}

			++sitesEqualToCoo;
			discardedTerms.push_back(i);

			OneOperatorSpecType opspec(siteSplit.root);
			const PsimagLite::String destKet = tmp + "*" + ket;
			OperatorType* op = new OperatorType(aux_.aoe(). model().
			                                    naturalOperator(opspec.label,
			                                                    0, // FIXME TODO SDHS Immm
			                                                    opspec.dof));
			if (opspec.transpose) op->transpose();

			oneOperator(destKet, ket, *op, site);
			ket = "|!m" + tmp + "*" + ket;
			factor_ = 1;
			strFactor_ = "";
			delete op;
			op = 0;
		}

		// discarded terms and inversion
		const SizeType nnew = newVstr.size();
		vStr_.clear();
		for (SizeType i = 0; i < nnew; ++i)
			vStr_.push_back(newVstr[nnew - i - 1]);
		vStr_[nnew - 1] = ket;

		finalized_ = true;
	}

	PsimagLite::String toString() const
	{
		PsimagLite::String s;
		const SizeType n = vStr_.size();
		if (n == 0)
			err("toString returns empty\n");

		PsimagLite::String f  = (strFactor_ != "") ? strFactor_ + "*" : "";

		for (SizeType i = 0; i < n - 1; ++i)
			s += vStr_[i] + "*";
		return f + s + vStr_[n - 1];
	}

	void setString(PsimagLite::String str)
	{
		if (vStr_.size() != 1)
			err("TermForTargetingExpression::setString\n");
		vStr_[0] = str;
	}

	SizeType size() const { return vStr_.size(); }

	bool finalized() const { return finalized_; }

	int pIndex() const
	{
		if (vStr_.size() != 1) return -1;
		const PsimagLite::String str = vStr_[0];
		SizeType last = str.length();
		if (last < 4) return -1;
		--last;
		if (str.substr(0, 2) == "|P" && str[last] == '>')
			return PsimagLite::atoi(str.substr(2, last - 2));
		return -1;
	}

private:

	bool siteCanBeApplied(SizeType site) const
	{
		const SizeType currentCoO = getCurrentCoO();
		const SizeType linSize = aux_.aoe().model().superGeometry().numberOfSites();
		const bool b1 = (site == 0 && currentCoO == 1 &&
		                 aux_.direction() == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);
		const bool b2 = (site == linSize - 1 && currentCoO == linSize - 2 &&
		                 aux_.direction() == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);

		return (b1 || b2 || site == currentCoO);
	}

	void oneOperator(PsimagLite::String destKet,
	                 PsimagLite::String srcKet,
	                 const OperatorType& op,
	                 SizeType site)
	{
		assert(siteCanBeApplied(site));
		const VectorWithOffsetType& srcVwo = aux_.getCurrentVectorConst(srcKet);
		PsimagLite::String internalName = aux_.createTemporaryVector(destKet);
		VectorWithOffsetType& destVwo = aux_.getCurrentVectorNonConst(internalName);
		applyInSitu(destVwo, srcVwo, site, op);
		destVwo *= factor_;
	}

	// returns A|src1>
	void applyInSitu(VectorWithOffsetType& dest,
	                 const VectorWithOffsetType& src1,
	                 SizeType site,
	                 const OperatorType& A)
	{
		const SizeType splitSize = aux_.aoe().model().hilbertSize(site);

		typename PsimagLite::Vector<bool>::Type oddElectrons;
		aux_.aoe().model().findOddElectronsOfOneSite(oddElectrons, site);
		FermionSign fs(aux_.lrs().left(), oddElectrons);
		bool b1 = (site == 0);
		SizeType n = aux_.aoe().model().superGeometry().numberOfSites();
		assert(n > 2);
		bool b2 = (site == n - 1);
		BorderEnumType border = (b1 || b2) ? BorderEnumType::BORDER_YES
		                                   : BorderEnumType::BORDER_NO;
		aux_.aoe().applyOpLocal()(dest, src1, A, fs, splitSize, aux_.direction(), border);
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

	void finalMultReal()
	{
		assert(PsimagLite::imag(factor_) == 0);

		const RealType f = PsimagLite::real(factor_);

		if (f < 0)
			strFactor_ = "(" + ttos(f) + ")";
		else
			strFactor_ = ttos(f);

		if (f == 1) strFactor_ = "";
	}

	void finalMultImag()
	{
		const RealType freal = PsimagLite::real(factor_);
		const RealType fimag = PsimagLite::imag(factor_);

		strFactor_ = "(" + ttos(freal) + "+" + ttos(fimag) + "i)";

		if (freal == 1 && fimag == 0) strFactor_ = "";
	}

	static bool isGlobalOperator(PsimagLite::String op)
	{
		return (op == "TimeEvolve@0@");
	}

	bool finalized_;
	const AuxiliaryType& aux_;
	ComplexOrRealType factor_;
	PsimagLite::String strFactor_;
	VectorStringType vStr_;
	NonLocalForTargetingExpressionType nonLocal_;
};
}
#endif // TERMFORTARGETINGEXPRESSION_H
