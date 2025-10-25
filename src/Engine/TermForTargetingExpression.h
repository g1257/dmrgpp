#ifndef TERMFORTARGETINGEXPRESSION_H
#define TERMFORTARGETINGEXPRESSION_H
#include "AuxForTargetingExpression.h"
#include "NonLocalForTargetingExpression.h"
#include "OneOperatorSpec.h"
#include "Vector.h"
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
	using KetType = KetForTargetingExpression<ComplexOrRealType>;

	TermForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false)
	    , aux_(aux)
	    , nonLocal_(aux)
	{
	}

	TermForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	    : finalized_(false)
	    , aux_(aux)
	    , nonLocal_(aux)
	{
		if (str.substr(0, 1) == "|") {
			ket_.set(str);
		} else {
			ops_.push_back(str);
		}
	}

	TermForTargetingExpression& operator=(const TermForTargetingExpression& other) = delete;

	TermForTargetingExpression(const TermForTargetingExpression& other)
	    : finalized_(other.finalized_)
	    , aux_(other.aux_)
	    , nonLocal_(other.aux_)
	{
		ops_ = other.ops_;
		ket_ = other.ket_;
	}

	void assignAndDestroy(TermForTargetingExpression& other)
	{
		finalized_ = other.finalized_;
		ops_ = other.ops_;
		ket_ = other.ket_;
	}

	void multiply(const TermForTargetingExpression& other)
	{
		const SizeType n = other.ops_.size();
		for (SizeType i = 0; i < n; ++i)
			ops_.push_back(other.ops_[i]);

		bool thisNoKet = this->ket_.name().empty();
		bool otherNoKet = other.ket().name().empty();

		if (!otherNoKet) {
			if (!thisNoKet) {
				err("Cannot multiply two vectors, only matrix times vector\n");
			}

			ket_ = other.ket();
		}

		if (n > 0) finalized_ = false;

	}

	void multiply(const ComplexOrRealType& val)
	{
		ket_.multiply(val);
	}

	void sum(const TermForTargetingExpression& other)
	{
		ket_.sum(other.ket());
	}

	void setFactor(const ComplexOrRealType& val)
	{
		ket_.setFactor(val);
	}

	void finalize()
	{
		if (finalized_)
			return;

		SizeType n = ops_.size();
		if (n == 0) {
			finalized_ = true;
			return;
		}

		SizeType sitesEqualToCoo = 0;
		VectorSizeType discardedTerms;

		VectorStringType newVstr;
		std::string reading_buffer;
		for (SizeType ii = 0; ii < n; ++ii) {
			const SizeType i = n - ii - 1; // read vector backwards
			PsimagLite::String tmp = ops_[i];
			reading_buffer += ops_[i] + "*";

			// it's a matrix
			if (ii != 0) {
				newVstr.push_back(tmp);
				continue; // apply in order only IMPORTANT
				// the last is the ket
			}

			SiteSplitType siteSplit = OneOperatorSpecType::extractSiteIfAny(tmp);
			if (isGlobalOperator(tmp)) {
				nonLocal_.timeEvolve(tmp, siteSplit, ket_.name(), aux_.currentCoO());
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

	PsimagLite::String toString() const
	{
		return ket_.toString(ops_);
	}

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

	const KetType& ket() const { return ket_; }

	bool finalized() const { return finalized_; }

	int pIndex() const
	{
		if (ops_.size() > 0)
			return -1;
		return ket_.pIndex();
	}

	bool isPureKet() const
	{
		return (ops_.empty());
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
	VectorStringType ops_;
	KetType ket_;
	NonLocalForTargetingExpressionType nonLocal_;
};
}
#endif // TERMFORTARGETINGEXPRESSION_H
