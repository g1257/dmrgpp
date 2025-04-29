#ifndef ALGEBRAFORTARGETINGEXPRESSION_H
#define ALGEBRAFORTARGETINGEXPRESSION_H

#include "AuxForTargetingExpression.h"
#include "GetBraOrKet.h"
#include "PackIndices.h"
#include "TermForTargetingExpression.h"
#include "Vector.h"

namespace Dmrg
{

template <typename TargetingBaseType>
class AlgebraForTargetingExpression
{

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

	AlgebraForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false)
	    , aux_(aux)
	{
	}

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

	AlgebraForTargetingExpression(const AlgebraForTargetingExpression&) = delete;

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

	bool isEmpty() const { return (terms_.size() == 0); }

	void finalize()
	{
		const SizeType nterms = terms_.size();
		for (SizeType i = 0; i < nterms; ++i)
			terms_[i]->finalize();
	}

	bool finalized() const { return finalized_; }

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

	int pIndex() const
	{
		if (terms_.size() != 1)
			return -1;
		return terms_[0]->pIndex();
	}

	bool hasSummationKetAndNoMult() const
	{
		const SizeType nterms = terms_.size();
		bool summation = false;
		bool mult = false;
		for (SizeType i = 0; i < nterms; ++i) {
			if (terms_[i]->size() != 1)
				continue;
			if (terms_[i]->toString().substr(0, 3) == "|!a")
				summation = true;
			if (terms_[i]->toString().substr(0, 3) == "|!m")
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
			if (!isTermSimplifiable(i)) {
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

				simpleTerms[j].factor += terms_[i]->factor();
				term_mapping[i] = -2; //ignore because i-th summed to j-th
				break;
			}

			if (term_mapping[i] == -1) {
				SimpleTerm st;
				st.index = i;
				st.enabled = true;
				st.factor = terms_[i]->factor();
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
		VectorTermType termsNew;
		const SizeType nterms = terms_.size();
		VectorBoolType indices(nterms, false);
		SizeType survivingTermIndex = nterms + 1000;
		SizeType survivingTerms = 0;

		for (SizeType i = 0; i < nterms; ++i) {
			if (!isTermSimplifiable(i)) {
				termsNew.push_back(terms_[i]);
				indices[i] = true;
				continue;
			}

			int pIndex = terms_[i]->pIndex();
			if (pIndex < 0) {
				err("Internal Error: pIndex < for for " +
				    terms_[i]->toString() + "\n");
			}

			if (survivingTerms == 0) {
				survivingTermIndex = termsNew.size();
				termsNew.push_back(terms_[i]);
				indices[i] = true;
				name = "|!aP" + ttos(pIndex);
			} else {
				name += "pP" + ttos(pIndex);
			}

			++survivingTerms;
		}

		if (survivingTerms < 2)
			return;

		name += ">";
		termsNew[survivingTermIndex]->setString(name);

		for (SizeType i = 0; i < nterms; ++i) {
			if (indices[i])
				continue;
			delete terms_[i];
			terms_[i] = nullptr;
		}

		terms_.swap(termsNew);
	}

	bool isTermSimplifiable(SizeType i) const
	{
		if (terms_[i]->size() != 1 || !terms_[i]->finalized())
			return false;
		PsimagLite::String str = terms_[i]->component(0);
		return (str.substr(0, 2) != "|!");
	}

	// important: needs operator=
	bool finalized_;
	const AuxiliaryType& aux_;
	VectorTermType terms_;
};

}
#endif // ALGEBRAFORTARGETINGEXPRESSION_H
