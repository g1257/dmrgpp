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

template<typename TargetingBaseType>
struct AuxForTargetingExpression {

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;

	AuxForTargetingExpression(const ApplyOperatorExpressionType& aoe_,
	                          const ModelType& model_,
	                          const LeftRightSuperType& lrs_,
	                          const VectorWithOffsetType& gs_,
	                          const VectorVectorWithOffsetType& pvectors_,
	                          ProgramGlobals::DirectionEnum dir)
	    : aoe(aoe_), model(model_), lrs(lrs_), gs(gs_), pvectors(pvectors_), direction(dir)
	{}

	const ApplyOperatorExpressionType& aoe;
	const ModelType& model;
	const LeftRightSuperType lrs;
	const VectorWithOffsetType& gs;
	const VectorVectorWithOffsetType& pvectors;
	ProgramGlobals::DirectionEnum direction;
};

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
	typedef OneOperatorSpec OneOperatorSpecType;
	typedef typename PsimagLite::Vector<OneOperatorSpecType*>::Type VectorOneOperatorSpecType;
	typedef PsimagLite::Vector<int>::Type VectorIntType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename OperatorType::StorageType SparseMatrixType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename ApplyOperatorExpressionType::BorderEnumType BorderEnumType;

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
		fullVector_ += otherCopy.fullVector_;
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
				*vwo = fullVector_;
				fullVector_.clear();
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
			fullVector_ = factor_*fullVector_;
		factor_ = 1.0;

		if (vwo) {
			*vwo = fullVector_;
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
		SizeType n = ops.size();
		assert(n == sites.size());
		OperatorType* op0 = 0;
		OperatorType* op1  =0;
		switch (n) {
		case 1:
			op0 = new OperatorType(aux_.model.naturalOperator(ops[0]->label,
			                       0, // FIXME TODO SDHS Immm
			                       ops[0]->dof));
			oneOperator(ket, *op0, sites[0]);
			delete op0;
			op0 = 0;
			break;

		case 2:
			op0 =  new OperatorType(aux_.model.naturalOperator(ops[0]->label,
			                        0, // FIXME TODO SDHS Immm
			                        ops[0]->dof));
			op1 = new OperatorType(aux_.model.naturalOperator(ops[1]->label,
			                       0, // FIXME TODO SDHS Immm
			                       ops[1]->dof));
			twoOperators(ket, *op0, sites[0], *op1, sites[1]);
			delete op0;
			delete op1;
			op0 = op1 = 0;
			break;

		default:
			err("finalizeInternal: Only 1 or 2 operators supported for now\n");
			break;
		}
	}

	void oneOperator(PsimagLite::String ket, const OperatorType& op, SizeType site)
	{
		SizeType currentCoO = getCurrentCoO();

		if (site == currentCoO) {
			const VectorWithOffsetType& srcVwo = getCurrentVector(ket);
			applyInSitu(srcVwo, site, op);
			return;
		}

		err("oneOperator unimplemented\n");
		// Fetch ket at coo site --> into vec
		// Apply op to vec ---> vec2
		// Bring vec2 back to current coo
	}

	void twoOperators(PsimagLite::String ket,
	                  const OperatorType& op1,
	                  SizeType site1,
	                  const OperatorType& op2,
	                  SizeType site2)
	{
		if (site1 == site2) {
			oneOperator(ket, op1*op2, site1);
			return;
		}

		SizeType currentCoO = getCurrentCoO();
		const VectorWithOffsetType& srcVwo = getCurrentVector(ket);
		if (site1 == currentCoO) {
			applyInSitu(srcVwo, site1, op1);
			oneOperator(ket, op2, site2);
			return;
		}

		if (site2 == currentCoO) {
			applyInSitu(srcVwo, site2, op2);
			oneOperator(ket, op1, site1);
			return;
		}

		err("twoOperators unimplemented\n");
		// neither is current CoO
		// Fetch ket at coo site1 --> into vec
		// Apply op1 to vec ---> vec2
		// Move to coo site2 --> vec3
		// Apply op2 to vec3 --> vec4
		// Bring vec4 back to current coo
	}

	// returns A|src1>
	void applyInSitu(const VectorWithOffsetType& src1,
	                 SizeType site,
	                 const OperatorType& A)
	{
		err("applyInSitu unimplemented\n");

		typename PsimagLite::Vector<bool>::Type oddElectrons;
		aux_.model.findOddElectronsOfOneSite(oddElectrons,site);
		FermionSign fs(aux_.lrs.left(), oddElectrons);
		bool b1 = (site == 1 && aux_.direction == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON);
		SizeType n = aux_.model.geometry().numberOfSites();
		assert(n > 2);
		bool b2 = (site == n - 2 && aux_.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);
		BorderEnumType border = (b1 || b2) ? BorderEnumType::BORDER_YES
		                                   : BorderEnumType::BORDER_NO;
		aux_.aoe.applyOpLocal()(fullVector_, src1, A, fs, aux_.direction, border);
	}

	const VectorWithOffsetType& getCurrentVector(PsimagLite::String braOrKet) const
	{
		GetBraOrKet getBraOrKet(braOrKet);

		SizeType ind = getBraOrKet();

		if (ind > 0 && ind - 1 >= aux_.pvectors.size())
			err("getVector: out of range for " + braOrKet + "\n");

		return (ind == 0) ? aux_.gs : aux_.pvectors[ind - 1];
	}

	SizeType getCurrentCoO() const
	{
		const LeftRightSuperType& lrs = aux_.lrs;
		const SizeType systemBlockSize = lrs.left().block().size();
		assert(systemBlockSize > 0);
		const int maxSystemSite = lrs.left().block()[systemBlockSize - 1];
		return (aux_.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? maxSystemSite
		                                                                        : maxSystemSite + 1;

	}

	bool finalized_;
	VectorStringType vStr_;
	VectorWithOffsetType fullVector_;
	ComplexOrRealType factor_;
	const AuxiliaryType& aux_;
};

template<typename TargetingBaseType>
class SpecForTargetingExpression {

public:

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef AlgebraForTargetingExpression<TargetingBaseType> AlgebraType;
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
