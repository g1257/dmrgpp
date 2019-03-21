#ifndef SPECFORTARGETINGEXPRESSION_H
#define SPECFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "OneOperatorSpec.h"
#include "CanonicalExpression.h"
#include <numeric>
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "PackIndices.h"

// All this isn't efficient!!

namespace Dmrg {

template<typename VectorWithOffsetType, typename ModelType>
struct AuxForTargetingExpression {

	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;

	AuxForTargetingExpression(const ModelType& model_,
	                          const LeftRightSuperType& lrs_,
	                          const VectorWithOffsetType& gs_,
	                          const VectorVectorWithOffsetType& pvectors_,
	                          ProgramGlobals::DirectionEnum dir)
	    : model(model_), lrs(lrs_), gs(gs_), pvectors(pvectors_), direction(dir)
	{}

	const ModelType& model;
	const LeftRightSuperType lrs;
	const VectorWithOffsetType& gs;
	const VectorVectorWithOffsetType& pvectors;
	ProgramGlobals::DirectionEnum direction;
};

template<typename VectorWithOffsetType, typename ModelType>
class AlgebraForTargetingExpression {

public:

	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef AuxForTargetingExpression<VectorWithOffsetType, ModelType> AuxiliaryType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename ModelType::ModelHelperType ModelHelperType;
	typedef typename ModelHelperType::LeftRightSuperType LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::OperatorsType OperatorsType;
	typedef typename OperatorsType::OperatorType OperatorType;
	typedef OneOperatorSpec OneOperatorSpecType;
	typedef typename PsimagLite::Vector<OneOperatorSpecType*>::Type VectorOneOperatorSpecType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;

	class GetOperator {

	public:

		GetOperator(SizeType index,
		            const BasisWithOperatorsType& basis,
		            bool transpose)
		    : index_(index), m_(0), owner_(false), fs_(1)
		{
			const OperatorType& op = basis.getOperatorByIndex(index);
			fs_ = (op.fermionOrBoson == ProgramGlobals::FermionOrBosonEnum::FERMION) ? -1 : 1;
			if (!transpose) {
				m_ = new SparseMatrixType();
				SparseMatrixType& mm = const_cast<SparseMatrixType&>(*m_);
				transposeConjugate(mm, op.data);
				owner_ = true;
			} else {
				m_ = &op.data;
			}
		}

		~GetOperator()
		{
			if (!owner_) return;
			delete m_;
			m_ = 0;
		}

		const SparseMatrixType& operator()() const
		{
			return *m_;
		}

		const int& fermionicSign() const { return fs_; }

	private:

		SizeType index_;
		SparseMatrixType const*  m_;
		bool owner_;
		int fs_;
	};

	AlgebraForTargetingExpression(const AuxiliaryType& aux)
	    : finalized_(false), factor_(1.0), aux_(aux) {}

	AlgebraForTargetingExpression(PsimagLite::String str, const AuxiliaryType& aux)
	    : finalized_(false),
	      vStr_(1, str),
	      factor_(1.0),
	      aux_(aux)
	{}

	// AlgebraForTargetingExpression(const AlgebraForTargetingExpression&) = delete;

	AlgebraForTargetingExpression& operator=(const AlgebraForTargetingExpression& other)
	{
		finalized_ = other.finalized_;
		vStr_ = other.vStr_;
		fullVector_ = other.fullVector_;
		factor_ = other.factor_;
		return *this;
	}

	AlgebraForTargetingExpression& operator+=(const AlgebraForTargetingExpression& other)
	{
		AlgebraForTargetingExpression otherCopy = other;
		otherCopy.finalize(0);
		finalize(0);
		fullVector_ += fullVector_;
		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const AlgebraForTargetingExpression& other)
	{
		if (other.finalized_ || finalized_)
			err("AlgebraForTargetingExpression: Two finalized terms cannot be multiplied\n");

		vStr_.insert(vStr_.end(), other.vStr_.begin(), other.vStr_.end());

		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const ComplexOrRealType& scalar)
	{
		factor_ *= scalar;
		return *this;
	}

	const VectorStringType& vStr() const { return vStr_; }

	PsimagLite::String toString() const
	{
		PsimagLite::String s;
		std::accumulate(vStr_.begin(), vStr_.end(), s);
		return s;
	}

	bool finalized() const { return finalized_; }

	void finalize(VectorWithOffsetType* vwo)
	{
		if (finalized_) {
			if (vwo) {
				vwo->fromFull(fullVector_, aux_.lrs.super());
				std::fill(fullVector_.begin(), fullVector_.end(), 0.0);
			}

			return;
		}

		SizeType n = vStr_.size();
		if (n == 0)
			err("AlgebraForTargetingExpression: Cannot finalize an empty object\n");

		SizeType opsSize = (n == 1) ? 1 : n - 1;
		VectorOneOperatorSpecType ops(opsSize, 0);
		VectorIntType sites(opsSize, -1);
		SizeType j = 0;
		PsimagLite::String ket;

		for (SizeType i = 0; i < n; ++i) {
			PsimagLite::String tmp = vStr_[i];
			if (tmp[0] == '|') { // it's a vector
				if (ket != "")
					err("More than one ket found in " + toString() + "\n");
				ket = tmp;
				if (i + 1 != n)
					err("Vector is not at the end in " + toString() + "\n");

				continue; // == break;
			}

			// it's a matrix
			assert(j < sites.size());
			sites[j] = OneOperatorSpecType::extractSiteIfAny(tmp);
			assert(j < ops.size());
			ops[j] = new OneOperatorSpecType(tmp);
			++j;
		}

		if (n > 1)
			finalizeInternal(ket, ops, sites);

		for (SizeType i = 0; i < opsSize; ++i) {
			delete ops[i];
			ops[i] = 0;
		}

		if (factor_ != 1.0)
			fullVector_ *= factor_;
		factor_ = 1.0;

		if (vwo) {
			vwo->fromFull(fullVector_, aux_.lrs.super());
			fullVector_.clear();
		}

		finalized_ = true;
		vStr_.clear();
	}

private:


	void finalizeInternal(PsimagLite::String ket,
	                      const VectorOneOperatorSpecType& ops,
	                      const VectorIntType& sites)
	{
//		SizeType n = ops.size();
//		assert(n == sites.size());
//		switch (n) {
//		case 1:

//			break;
//		default:
//			break;
//		}
	}

	void checkSites(const VectorIntType& sites) const
	{
		SizeType n = sites.size();
		for (SizeType i = 1; i < n; ++i) {
			if (sites[i] >= sites[i - 1]) continue;
			err(PsimagLite::String("SpecForTargetingExpression: Sites must be ") +
			    " ordered increasingly in expression " + toString());
		}
	}

	bool finalized_;
	VectorStringType vStr_;
	VectorType fullVector_;
	ComplexOrRealType factor_;
	const AuxiliaryType& aux_;
};

template<typename VectorWithOffsetType, typename ModelType>
class SpecForTargetingExpression {

public:

	typedef AlgebraForTargetingExpression<VectorWithOffsetType, ModelType> AlgebraType;
	typedef AlgebraType ResultType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename AlgebraType::AuxiliaryType AuxiliaryType;

	AlgebraType operator()(PsimagLite::String str, const AuxiliaryType& aux) const
	{
		return AlgebraType(str, aux);
	}

	static bool metaEqual(const ResultType&, const ResultType&) { return true; }

	static bool isEmpty(const ResultType& term)
	{
		if (term.finalized()) return false;
		return (term.vStr().size() == 0);
	}
};
}
#endif // SPECFORTARGETINGEXPRESSION_H
