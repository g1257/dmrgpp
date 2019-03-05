#ifndef SPECFORTARGETINGEXPRESSION_H
#define SPECFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "OneOperatorSpec.h"
#include "CanonicalExpression.h"
#include <numeric>
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "PackIndices.h"

namespace Dmrg {

template<typename VectorWithOffsetType, typename ModelType>
struct AuxForTargetingExpression {

	typedef typename PsimagLite::Vector<VectorWithOffsetType>::Type VectorVectorWithOffsetType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;

	AuxForTargetingExpression(const ModelType& model_,
	                          const LeftRightSuperType& lrs_,
	                          const VectorWithOffsetType& gs_,
	                          const VectorVectorWithOffsetType& pvectors_)
	    : model(model_), lrs(lrs_), gs(gs_), pvectors(pvectors_)
	{}

	const ModelType& model;
	const LeftRightSuperType lrs;
	const VectorWithOffsetType& gs;
	const VectorVectorWithOffsetType& pvectors;
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
	typedef typename VectorWithOffsetType::VectorType VectorType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename OperatorType::StorageType SparseMatrixType;

	const VectorWithOffsetType& toVectorWithOffset() const { return data_; }

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
		data_ = other.data_;
		factor_ = other.factor_;
		return *this;
	}

	AlgebraForTargetingExpression& operator+=(const AlgebraForTargetingExpression& other)
	{
		AlgebraForTargetingExpression otherCopy = other;
		otherCopy.finalize();
		finalize();
		data_ += other.toVectorWithOffset();
		return *this;
	}

	AlgebraForTargetingExpression& operator*=(const AlgebraForTargetingExpression& other)
	{
		if (other.finalized_ || finalized_)
			err("AlgebraForTargetingExpression: Two finalized terms cannot be multiplied\n");

		vStr_.insert(vStr_.begin(), other.vStr_.begin(), other.vStr_.end());

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

	void finalize()
	{
		if (finalized_) return;

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

		if (data_.size() != 0)
			err("AlgebraForTargetingExpression: data already set !?\n");

		data_ = getVector(ket);
		data_ = factor_*data_;
		factor_ = 1.0;

		if (n > 1)
			finalizeInternal(ops, sites);

		for (SizeType i = 0; i < opsSize; ++i) {
			delete ops[i];
			ops[i] = 0;
		}

		finalized_ = true;
		vStr_.clear();
	}

private:

	void finalizeInternal(const VectorOneOperatorSpecType& ops,
	                      const VectorIntType& sites)
	{
		SizeType sectors = data_.sectors();
		for (SizeType i = 0; i < sectors; ++i) {
			finalizeInternal(ops, sites, data_.sector(i));
		}
	}

	void finalizeInternal(const VectorOneOperatorSpecType& ops,
	                      const VectorIntType& sites,
	                      SizeType iSector)
	{
		SizeType n = ops.size();
		const SizeType opsPerSite = aux_.model.modelLinks().cm().size();
		const LeftRightSuperType& lrs = aux_.lrs;
		const SizeType systemBlockSize = lrs.left().block().size();
		assert(systemBlockSize > 0);
		const int maxSystemSite = lrs.left().block()[systemBlockSize - 1];
		VectorType v(lrs.super().size(), 0.0);

		for (SizeType i = 0; i < n; ++i) {
			SizeType j = n - i - 1;
			SizeType index = aux_.model.modelLinks().nameDofToIndex(ops[j]->label,
			                                                        ops[j]->dof);
			assert(j < sites.size());

			if (sites[j] <= maxSystemSite) { // in system
				index += opsPerSite*sites[j];

				const SparseMatrixType& m = lrs.left().getOperatorByIndex(index).data;
				multiplySys(v, iSector, m, ops[j]->transpose);
			} else { // in environ
				SizeType siteReverse = aux_.model.geometry().numberOfSites() - sites[j] - 1;
				index += opsPerSite*siteReverse;
				const SparseMatrixType& m = lrs.right().getOperatorByIndex(index).data;
				multiplyEnv(v, iSector, m, ops[j]->transpose);
			}

			data_.fromFull(v, aux_.lrs.super()); // doesn't work in general FIXME
			// convert data_ from the beginning into BlockOffDiag.h FIXME TODO
			err("AlgebraForTargetingExpression::finalizeInternal() not implemented\n");

			std::fill(v.begin(), v.end(), 0.0);
		}
	}

	void multiplySys(VectorType& v,
	                 SizeType iSector,
	                 const SparseMatrixType& m,
	                 bool transpose)
	{
		PackIndicesType pack(aux_.lrs.left().size());
		SizeType offset = data_.offset(iSector);
		SizeType total = data_.effectiveSize(iSector);
		for (SizeType i = 0; i < total; ++i) {
			SizeType s = 0;
			SizeType e = 0;
			pack.unpack(s, e, aux_.lrs.super().permutation(i + offset));
			const SizeType start = m.getRowPtr(s);
			const SizeType end = m.getRowPtr(s + 1);
			for (SizeType k = start; k < end; ++k) {
				const SizeType sprime = m.getCol(k);
				const SizeType j = pack.pack(sprime, e, aux_.lrs.super().permutationInverse());
				assert(j < v.size());
				v[j] += m.getValue(k)*data_.fastAccess(iSector, i);
			}
		}
	}

	void multiplyEnv(VectorType& v,
	                 SizeType iSector,
	                 const SparseMatrixType& m,
	                 bool transpose)
	{
		err("AlgebraForTargetingExpression::multiplyEnv() not implemented\n");
	}

	const VectorWithOffsetType& getVector(PsimagLite::String braOrKet) const
	{
		GetBraOrKet getBraOrKet(braOrKet);

		SizeType ind = getBraOrKet();

		if (ind > 0 && ind - 1 >= aux_.pvectors.size())
			err("getVector: out of range for " + braOrKet + "\n");

		return (ind == 0) ? aux_.gs : aux_.pvectors[ind - 1];
	}

	bool finalized_;
	VectorStringType vStr_;
	VectorWithOffsetType data_;
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
