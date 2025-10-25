#ifndef TERMFORTARGETINGEXPRESSION_H
#define TERMFORTARGETINGEXPRESSION_H
#include "AuxForTargetingExpression.h"
#include "NonLocalForTargetingExpression.h"
#include "OneOperatorSpec.h"
#include "Vector.h"
#include "FactorForTargetingExpression.hh"
#include "KetForTargetingExpression.hh"

namespace Dmrg
{

template <typename TargetingBaseType>
class TermForTargetingExpression
{

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
	using VectorType = std::vector<ComplexOrRealType>;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef AuxForTargetingExpression<TargetingBaseType> AuxiliaryType;
	typedef typename AuxiliaryType::PvectorsType PvectorsType;
	typedef PsimagLite::OneOperatorSpec OneOperatorSpecType;
	typedef typename PsimagLite::Vector<OneOperatorSpecType*>::Type VectorOneOperatorSpecType;
	typedef typename OneOperatorSpecType::SiteSplit SiteSplitType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename ApplyOperatorExpressionType::BorderEnumType BorderEnumType;
	typedef NonLocalForTargetingExpression<TargetingBaseType> NonLocalForTargetingExpressionType;
	using FactorForTargetingExpressionType = FactorForTargetingExpression<ComplexOrRealType>;

	TermForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false)
	    , aux_(aux)
	    , factor_(1.0)
	    , nonLocal_(aux)
	{
	}

	// TermForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	//     : finalized_(false)
	//     , aux_(aux)
	//     , type_(E)
	//     , factor_(1.0)
	//     , vStr_(1, str)
	//     , nonLocal_(aux)
	// {
	// }

	TermForTargetingExpression& operator=(const TermForTargetingExpression& other) = delete;

	TermForTargetingExpression(const TermForTargetingExpression& other) = delete;

	void assignAndDestroy(TermForTargetingExpression& other)
	{
		finalized_ = other.finalized_;
		factor_ = other.factor_;
		ops_ = other.ops_;
		ket_ = other.ket_;
	}

	void multiply(const TermForTargetingExpression& other)
	{
		const SizeType n = other.ops_.size();
		for (SizeType i = 0; i < n; ++i)
			ops_.push_back(other.ops_[i]);
	}

	void multiply(const ComplexOrRealType& val)
	{
		factor_.multiply(val);
	}

	void sum(const TermForTargetingExpression& other)
	{
		err("sum of terms unimplemented\n");
	}

	void setFactor(const ComplexOrRealType& val)
	{
		factor_.set(val);
	}

	void setFactor(const VectorType& vec)
	{
		factor_.set(vec);
	}

	void finalize()
	{
		if (finalized_)
			return;

		SizeType n = ops_.size();
		if (n == 0)
			err("AlgebraForTargetingExpression: Cannot finalize an empty object\n");

		SizeType sitesEqualToCoo = 0;
		VectorSizeType discardedTerms;

		VectorStringType newVstr;
		std::string reading_buffer;
		for (SizeType ii = 0; ii < n; ++ii) {
			const SizeType i = n - ii - 1; // read vector backwards
			PsimagLite::String tmp = ops_[i];
			reading_buffer += ops_[i] + "*";
			if (tmp[0] == '|') { // it's a vector
				if (!ket_.name().empty())
					err("More than one ket found in " + reading_buffer + "\n");
				ket_.set(tmp);
				if (i + 1 != n)
					err("Vector is not at the end in " + reading_buffer + "\n");
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
				nonLocal_.timeEvolve(tmp, siteSplit, ket, aux_.currentCoO());
				newVstr.push_back(tmp);
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
			std::string ket_src = ket_.name();
			ket_.multiply(tmp);
			std::string ket_dest = ket_.name();
			OperatorType* op = new OperatorType(aux_.pVectors().aoe().model().naturalOperator(opspec.label,
			    0, // FIXME TODO SDHS Immm
			    opspec.dof));
			if (opspec.transpose)
				op->transpose();

			oneOperator(ket_dest, ket_src, *op, site);

			delete op;
			op = 0;
		}

		// discarded terms and inversion
		const SizeType nnew = newVstr.size();
		ops_.clear();
		for (SizeType i = 0; i < nnew; ++i)
			ops_.push_back(newVstr[nnew - i - 1]);

		finalized_ = true;
	}

	void setKet(const std::string& ket)
	{
		ket_.set(ket);
	}

	// Constant functions below

	// PsimagLite::String toString() const
	// {
	// 	PsimagLite::String s;
	// 	const SizeType n = ops_.size();
	// 	if (n == 0)
	// 		err("toString returns empty\n");

	// 	PsimagLite::String f = factor_.toString();

	// 	for (SizeType i = 0; i < n - 1; ++i)
	// 		s += ops_[i] + "*";
	// 	return f + s + ket_.toString();
	// }

	// void setString(PsimagLite::String str)
	// {
	// 	if (vStr_.size() != 1)
	// 		err("TermForTargetingExpression::setString\n");
	// 	vStr_[0] = str;
	// }

	// SizeType size() const { return ops_.size() + 1; }

	// const std::string& component(SizeType ind) const
	// {
	// 	assert(ind < vStr_.size());
	// 	return vStr_[ind];
	// }

	bool isSummable() const
	{
		if (!ops_.empty() || !finalized_)
			return false;

		return ket_.isSummable();
	}

	ComplexOrRealType factor() const { return factor_.value(); }

	const KetForTargetingExpression& ket() const { return ket_; }

	bool finalized() const { return finalized_; }

	int pIndex() const
	{
		if (ops_.size() > 0)
			return -1;
		return ket_.pIndex();
	}

private:

	static bool isGlobalOperator(PsimagLite::String opName)
	{
		return PvectorsType::PvectorType::isTimeEvolution(opName);
	}

	bool siteCanBeApplied(SizeType site) const
	{
		return (site == aux_.currentCoO());
	}

	void oneOperator(PsimagLite::String destKet,
	    PsimagLite::String srcKet,
	    const OperatorType& op,
	    SizeType site)
	{
		assert(siteCanBeApplied(site));
		const VectorWithOffsetType& srcVwo = aux_.pVectors().getCurrentVectorConst(srcKet);
		PsimagLite::String internalName = aux_.createTemporaryVector(destKet);
		VectorWithOffsetType& destVwo = aux_.getCurrentVectorNonConst(internalName);
		applyInSitu(destVwo, srcVwo, site, op);
	}

	// returns A|src1>
	void applyInSitu(VectorWithOffsetType& dest,
	    const VectorWithOffsetType& src1,
	    SizeType site,
	    const OperatorType& A)
	{
		const SizeType splitSize = aux_.pVectors().aoe().model().hilbertSize(site);

		typename PsimagLite::Vector<bool>::Type oddElectrons;
		aux_.pVectors().aoe().model().findOddElectronsOfOneSite(oddElectrons, site);
		FermionSign fs(aux_.pVectors().lrs().left(), oddElectrons);
		bool b1 = (site == 0);
		SizeType n = aux_.pVectors().aoe().model().superGeometry().numberOfSites();
		assert(n > 2);
		bool b2 = (site == n - 1);
		BorderEnumType border = (b1 || b2) ? BorderEnumType::BORDER_YES
						   : BorderEnumType::BORDER_NO;
		aux_.pVectors().aoe().applyOpLocal()(dest,
		    src1,
		    A,
		    fs,
		    splitSize,
		    aux_.direction(),
		    border);
	}

	bool finalized_;
	const AuxiliaryType& aux_;
	FactorForTargetingExpressionType factor_;
	VectorStringType ops_;
	KetForTargetingExpression ket_;
	NonLocalForTargetingExpressionType nonLocal_;
};
}
#endif // TERMFORTARGETINGEXPRESSION_H
