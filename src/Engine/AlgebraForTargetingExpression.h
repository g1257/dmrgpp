#ifndef ALGEBRAFORTARGETINGEXPRESSION_H
#define ALGEBRAFORTARGETINGEXPRESSION_H

#include "Vector.h"
#include "AuxForTargetingExpression.h"
#include "TermForTargetingExpression.h"
#include "GetBraOrKet.h"

namespace Dmrg {

template<typename TargetingBaseType>
class AlgebraForTargetingExpression {

public:

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef AuxForTargetingExpression<TargetingBaseType> AuxiliaryType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef TermForTargetingExpression<TargetingBaseType> TermType;
	typedef typename PsimagLite::Vector<TermType>::Type VectorTermType;

	AlgebraForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false), aux_(aux) {}

	AlgebraForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	    : finalized_(false), aux_(aux)
	{
		terms_.push_back(TermType(str, aux));
	}

	// AlgebraForTargetingExpression(const AlgebraForTargetingExpression&) = delete;

	AlgebraForTargetingExpression& operator=(const AlgebraForTargetingExpression& other)
	{
		finalized_ = other.finalized_;
		terms_ = other.terms_;
		return *this;
	}

	AlgebraForTargetingExpression& operator+=(const AlgebraForTargetingExpression& other)
	{
		if (other.terms_.size() != 1)
			err("Add only one term at a time\n");
		AlgebraForTargetingExpression otherCopy = other;
		// finalize each term of this and ...
		for (SizeType i = 0; i < terms_.size(); ++i)
			terms_[i].finalize();

		// ... and of otherCopy
		otherCopy.terms_[0].finalize();

		// add those of otherCopy to this and ...
		terms_.push_back(otherCopy.terms_[0]);

		// ... simplify if possible
		simplifyTerms();

		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const AlgebraForTargetingExpression& other)
	{
		if (terms_.size() != 1 && other.terms_.size() != 1)
			err("Only canonical expressions supported\n");

		terms_[0].multiply(other.terms_[0]);

		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const ComplexOrRealType& scalar)
	{
		if (terms_.size() != 1)
			err("scalar multiplication\n");

		terms_[0].multiply(scalar);

		return *this;
	}

	bool isEmpty() const
	{
		return (terms_.size() == 0);
	}

	void finalize()
	{
		const SizeType nterms = terms_.size();
		for (SizeType i = 0; i < nterms; ++i)
			terms_[i].finalize();
	}

	bool finalized() const { return finalized_; }

	PsimagLite::String toString() const
	{
		PsimagLite::String s;
		const SizeType n = terms_.size();
		if (n == 0)
			err("toString returns empty\n");

		for (SizeType i = 0; i < n - 1; ++i)
			s += terms_[i].toString() + "+";
		return s + terms_[n - 1].toString();
	}

	SizeType size() const { return terms_.size(); }

	int pIndex() const
	{
		if (terms_.size() != 1) return -1;
		return terms_[0].pIndex();
	}

	bool hasSummationKetAndNoMult() const
	{
		const SizeType nterms = terms_.size();
		bool summation = false;
		bool mult = false;
		for (SizeType i = 0; i < nterms; ++i) {
			if (terms_[i].size() != 1) continue;
			if (terms_[i].toString().substr(0, 3) == "|!a")
				summation = true;
			if (terms_[i].toString().substr(0, 3) == "|!m")
				mult = true;
		}

		return (summation && !mult);
	}

private:

	void simplifyTerms()
	{
		PsimagLite::String name;
		VectorTermType termsNew;
		const SizeType nterms = terms_.size();
		SizeType survivingTermIndex = nterms + 1000;
		SizeType survivingTerms = 0;
		for (SizeType i = 0; i < nterms; ++i) {
			if (!isTermSimplifiable(i)) {
				termsNew.push_back(terms_[i]);
				continue;
			}

			PsimagLite::GetBraOrKet ket(terms_[i].toString());
			const SizeType pIndex = ket.pIndex();

			if (survivingTerms == 0) {
				survivingTermIndex = termsNew.size();
				termsNew.push_back(terms_[i]);
				name = "|!aP" + ttos(pIndex);
			} else {
				name += "pP" + ttos(pIndex);
			}

			++survivingTerms;
		}

		if (survivingTerms < 2) return;

		name += ">";
		termsNew[survivingTermIndex].setString(name);
		terms_ = termsNew;
	}

	bool isTermSimplifiable(SizeType i) const
	{
		if (terms_[i].size() != 1 || !terms_[i].finalized())
			return false;
		PsimagLite::String str = terms_[i].toString();
		return (str.substr(0, 2) != "|!");
	}

	// important: needs operator=
	bool finalized_;
	const AuxiliaryType& aux_;
	VectorTermType terms_;
};

}
#endif // ALGEBRAFORTARGETINGEXPRESSION_H
