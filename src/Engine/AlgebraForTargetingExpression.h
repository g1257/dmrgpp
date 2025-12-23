#ifndef ALGEBRAFORTARGETINGEXPRESSION_H
#define ALGEBRAFORTARGETINGEXPRESSION_H

#include "AuxForTargetingExpression.h"
#include "GetBraOrKet.h"
#include "PackIndices.h"
#include "TermForTargetingExpression.h"
#include "Vector.h"

namespace Dmrg {

template <typename TargetingBaseType> class AlgebraForTargetingExpression {

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
	typedef typename PsimagLite::Vector<TermType*>::Type VectorTermType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	using KetType = typename TermType::KetType;

	AlgebraForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false)
	    , aux_(aux)
	{ }

	AlgebraForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	    : finalized_(false)
	    , aux_(aux)
	{
		terms_.push_back(new TermType(str, aux));
	}

	~AlgebraForTargetingExpression()
	{
		const SizeType n = terms_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete terms_[i];
			terms_[i] = nullptr;
		}
	}

	AlgebraForTargetingExpression(const AlgebraForTargetingExpression& other)
	    : finalized_(other.finalized_)
	    , aux_(other.aux_)
	{
		SizeType n = other.size();
		terms_.resize(n, nullptr);
		for (SizeType i = 0; i < n; ++i) {
			terms_[i] = new TermType(other.term(i));
		}
	}

	AlgebraForTargetingExpression& operator=(const AlgebraForTargetingExpression&) = delete;

	void assignAndDestroy(AlgebraForTargetingExpression& other)
	{
		if (terms_.size() != 0)
			err("assignAndDestroy: term isn't empty\n");

		finalized_ = other.finalized_;

		const SizeType n = other.terms_.size();

		terms_.resize(n, new TermType(aux_));

		for (SizeType i = 0; i < n; ++i)
			terms_[i]->assignAndDestroy(*(other.terms_[i]));
	}

	// other gets destroyed here
	void plus(AlgebraForTargetingExpression& other)
	{
		if (other.terms_.size() != 1)
			err("Add only one term at a time\n");
		// finalize each term of this and ...
		for (SizeType i = 0; i < terms_.size(); ++i)
			terms_[i]->finalize();

		// ... and of otherCopy
		other.terms_[0]->finalize();

		// add those of otherCopy to this and ...
		const SizeType ind = terms_.size();
		terms_.push_back(new TermType(aux_));
		terms_[ind]->assignAndDestroy(*(other.terms_[0]));

		// ... simplify if possible
		simplifyTerms();
	}

	void multiply(const AlgebraForTargetingExpression& other)
	{
		if (terms_.size() != 1 && other.terms_.size() != 1)
			err("Only canonical expressions supported\n");

		terms_[0]->multiply(*other.terms_[0]);
	}

	void multiplyScalar(const ComplexOrRealType& scalar)
	{
		if (terms_.size() != 1)
			err("scalar multiplication\n");

		terms_[0]->multiply(scalar);
	}

	void finalize()
	{
		const SizeType nterms = terms_.size();
		for (SizeType i = 0; i < nterms; ++i)
			terms_[i]->finalize();
	}

	bool finalized() const { return finalized_; }

	void setKet(SizeType ind, const std::string& ket)
	{
		assert(ind < terms_.size());
		terms_[ind]->setKet(ket);
	}

	void setFactor(SizeType ind, const ComplexOrRealType& val)
	{
		assert(ind < terms_.size());
		terms_[ind]->setFactor(val);
	}

	// Constant functions below

	bool isEmpty() const { return (terms_.size() == 0); }

	PsimagLite::String toString() const
	{
		PsimagLite::String s;
		const SizeType n = terms_.size();
		if (n == 0)
			err("toString returns empty\n");

		for (SizeType i = 0; i < n - 1; ++i)
			s += terms_[i]->toString() + "+";
		return s + terms_[n - 1]->toString();
	}

	SizeType size() const { return terms_.size(); }

	AuxiliaryType& aux()
	{
		AuxiliaryType* auxPtr = const_cast<AuxiliaryType*>(&aux_);
		return *auxPtr;
	}

	const TermType& term(SizeType ind) const
	{
		if (ind >= terms_.size()) {
			err("pIndex(): ind >= terms.size()\n");
		}

		if (!terms_[ind]) {
			err("Term is null\n");
		}

		return *terms_[ind];
	}

	bool hasSummationKetAndNoMult() const
	{
		const SizeType nterms = terms_.size();
		bool summation = false;
		bool mult = false;
		for (SizeType i = 0; i < nterms; ++i) {
			if (!terms_[i]->isPureKet())
				continue;
			const KetType& ket = terms_[i]->ket();
			if (ket.kind() == KetType::Kind::S)
				summation = true;
			if (ket.kind() == KetType::Kind::M)
				mult = true;
		}

		return (summation && !mult);
	}

private:

	void simplifyTerms()
	{
		simpleSimplification();
		complicatedSimplification();
	}

	void simpleSimplification()
	{
		struct SimpleTerm {
			SizeType index = 0;
			SizeType pIndex = 0;
			bool enabled = false;
			ComplexOrRealType factor = 0.;
		};

		const SizeType nterms = terms_.size();

		std::vector<SimpleTerm> simpleTerms;
		std::vector<int> term_mapping(nterms, -1);

		for (SizeType i = 0; i < nterms; ++i) {
			const TermType& term = *terms_[i];
			if (!term.isSummable()) {
				continue;
			}

			int pIndexInt = terms_[i]->pIndex();
			if (pIndexInt < 0) {
				continue;
			}

			SizeType pIndex = pIndexInt;

			for (SizeType j = 0; j < simpleTerms.size(); ++j) {
				if (!simpleTerms[j].enabled || simpleTerms[j].pIndex != pIndex) {
					continue;
				}

				simpleTerms[j].factor += terms_[i]->ket().factor();
				term_mapping[i] = -2; // ignore because i-th summed to j-th
				break;
			}

			if (term_mapping[i] == -1) {
				SimpleTerm st;
				st.index = i;
				st.enabled = true;
				st.factor = terms_[i]->ket().factor();
				st.pIndex = pIndex;
				term_mapping[i] = simpleTerms.size();
				simpleTerms.push_back(st);
			}
		}

		VectorTermType newterms;

		for (SizeType i = 0; i < nterms; ++i) {
			if (term_mapping[i] == -2) {
				delete terms_[i];
				terms_[i] = nullptr;
				continue;
			} else if (term_mapping[i] == -1) {
				newterms.push_back(terms_[i]);
			} else {
				assert(term_mapping[i] >= 0);
				SizeType ind = term_mapping[i];
				SizeType index = newterms.size();
				newterms.push_back(terms_[i]);
				assert(static_cast<SizeType>(term_mapping[i]) < simpleTerms.size());
				const SimpleTerm& st = simpleTerms[ind];
				assert(st.enabled);
				newterms[index]->setFactor(st.factor);
			}
		}

		terms_ = newterms;
	}

	void complicatedSimplification()
	{
		PsimagLite::String name;
		const SizeType nterms = terms_.size();
		VectorBoolType indices(nterms, false);
		TermType* combined_term = new TermType(aux_);
		// VectorType factors;

		for (SizeType i = 0; i < nterms; ++i) {
			const TermType& term = *terms_[i];
			if (!term.isSummable()) {
				indices[i] = true;
				continue;
			}

			combined_term->sum(term);
		}

		if (!combined_term->ket().canSumBeFinished()) {
			delete combined_term;
			combined_term = nullptr;
			return;
		}

		VectorTermType termsNew;
		for (SizeType i = 0; i < nterms; ++i) {
			if (indices[i]) {
				termsNew.push_back(terms_[i]);
			} else {
				delete terms_[i];
				terms_[i] = nullptr;
			}
		}

		termsNew.push_back(combined_term);

		terms_.swap(termsNew);
	}

	// important: needs operator=
	bool finalized_;
	const AuxiliaryType& aux_;
	VectorTermType terms_;
};

}
#endif // ALGEBRAFORTARGETINGEXPRESSION_H
