#ifndef LAZYALGEBRA_H
#define LAZYALGEBRA_H
#include "Vector.h"

namespace Dmrg {

template<typename OperatorType>
class LazyAlgebraFactor {

	typedef typename PsimagLite::Vector<OperatorType>::Type VectorOperatorType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename OperatorType::value_type ComplexOrRealType;

public:

	LazyAlgebraFactor() : ops_(1, OperatorType()), indices_(1, 1), overallFactor_(1.0)
	{}

	LazyAlgebraFactor(const OperatorType& op)
	    : ops_(1, op), indices_(1, 1), overallFactor_(1.0)
	{}

	const LazyAlgebraFactor& operator*=(const ComplexOrRealType& f)
	{
		overallFactor_ *= f;
		return *this;
	}

	SizeType metaDiff(const LazyAlgebraFactor& other) const
	{
		const SizeType n = indices_.size();
		if (other.indices_.size() != n) return 1;
		SizeType sum = 0;
		for (SizeType i = 0; i < n; ++i) {
			SizeType ind = indices_[i];
			if (ind != other.indices_[i]) return 2;
			if (ind == 0) continue;
			--ind;
			assert(ind < ops_.size() && ind < other.ops_.size());
			sum += ops_[ind].metaDiff(other.ops_[ind]);
		}

		return sum;
	}

	bool isEmpty() const
	{
		const SizeType n = indices_.size();
		for (SizeType i = 0; i < n; ++i) {
			SizeType ind = indices_[i];
			if (ind == 0) return false;
			--ind;
			assert(ind < ops_.size());
			if (ops_[ind].isEmpty()) continue;
		}

		return true;
	}

	friend LazyAlgebraFactor operator*(const LazyAlgebraFactor& a,
	                                   const LazyAlgebraFactor& b)
	{
		LazyAlgebraFactor c = a;
		c.overallFactor_ *= b.overallFactor_;
		const SizeType n = b.indices_.size();
		SizeType offset = c.ops_.size();
		for (SizeType i = 0; i < n; ++i) {
			SizeType ind = b.indices_[i];
			if (ind == 0) {
				c.indices_.push_back(0);
				continue;
			}

			--ind;
			assert(ind < b.ops_.size());
			c.ops_.push_back(b.ops_[ind]);
			c.indices_.push_back(ind + 1 + offset);
		}

		return c;
	}

private:

	VectorOperatorType ops_;
	VectorSizeType indices_;
	ComplexOrRealType overallFactor_;
};

template<typename OperatorType>
class LazyAlgebra {

public:

	typedef LazyAlgebraFactor<OperatorType> LazyAlgebraFactorType;
	typedef typename PsimagLite::Vector<LazyAlgebraFactorType>::Type VectorLazyAlgebraFactorType;
	typedef typename OperatorType::value_type ComplexOrRealType;

	LazyAlgebra() : factors_(1, OperatorType())
	{}

	LazyAlgebra(const OperatorType& op) : factors_(1, op)
	{}

	const LazyAlgebra& operator+=(const LazyAlgebra& f)
	{
		const SizeType n = f.factors_.size();
		for (SizeType i = 0; i < n; ++i)
			factors_.push_back(f.factors_[i]);
		return *this;
	}

	const LazyAlgebra& operator*=(const ComplexOrRealType& f)
	{
		const SizeType n = factors_.size();
		for (SizeType i = 0; i < n; ++i)
			factors_[i] *= f;
		return *this;
	}

	const LazyAlgebra& operator*=(const LazyAlgebra& other)
	{
		const SizeType n = factors_.size();
		const SizeType m = other.factors_.size();

		VectorLazyAlgebraFactorType newFactors;
		for (SizeType i = 0; i < n; ++i)
			for (SizeType j = 0; j < m; ++j)
					newFactors.push_back(factors_[i]*other.factors_[j]);

		factors_ = newFactors;
		return *this;
	}

	SizeType metaDiff(const LazyAlgebra& other) const
	{
		const SizeType n = factors_.size();
		if (other.factors_.size() != n) return 1;
		SizeType sum = 0;
		for (SizeType i = 0; i < n; ++i)
			sum += factors_[i].metaDiff(other.factors_[i]);

		return sum;
	}

	bool isEmpty() const
	{
		const SizeType n = factors_.size();
		for (SizeType i = 0; i < n; ++i)
			if (factors_[i].isEmpty()) continue;

		return true;
	}

private:

	VectorLazyAlgebraFactorType factors_;
};
}
#endif // LAZYALGEBRA_H
