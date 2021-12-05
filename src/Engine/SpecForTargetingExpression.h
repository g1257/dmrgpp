#ifndef SPECFORTARGETINGEXPRESSION_H
#define SPECFORTARGETINGEXPRESSION_H
#include "Vector.h"
#include "OneOperatorSpec.h"
#include "CanonicalExpression.h"
#include <numeric>
#include "GetBraOrKet.h"
#include "ProgramGlobals.h"
#include "PackIndices.h"
#include "AlgebraForTargetingExpression.h"

namespace Dmrg {

template<typename TargetingBaseType>
class SpecForTargetingExpression {

public:

	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef AlgebraForTargetingExpression<TargetingBaseType> AlgebraType;
	typedef AlgebraType ResultType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename AlgebraType::AuxiliaryType AuxiliaryType;
	typedef std::pair<PsimagLite::String, AuxiliaryType> PairStringAuxType;

	class AssignAndDestroy {

	public:

		AssignAndDestroy(const AlgebraType& empty)
		    : isValid_(true)
		{
			AlgebraType* emptyPtr = const_cast<AlgebraType*>(&empty);

			t_ = new AlgebraType(emptyPtr->aux());

			if (!empty.isEmpty())
				err("AssignAndDestroy: Cannot start with non empty algebra\n");
		}

		AssignAndDestroy(const PairStringAuxType& pair)
		    : t_(new AlgebraType(pair.first, pair.second)), isValid_(true)
		{}

		~AssignAndDestroy()
		{
			isValid_ = false;
			delete t_;
			t_ = nullptr;
		}

		void assignBackward(AlgebraType& t2)
		{
			makeSureItsValid();
			// t2 = t_;
			t2.assignAndDestroy(*t_);
			isValid_ = false;
		}

		void assign(AssignAndDestroy& t2)
		{
			t_->assignAndDestroy(*t2.t_); // t_ = t2
			t2.isValid_ = false;
		}

		void plusBackward(AlgebraType& t2)
		{
			makeSureItsValid();
			// t2 += t_;
			t2.plus(*t_);
			isValid_ = false;
		}

		void multiply(AssignAndDestroy& t2) const
		{
			makeSureItsValid();
			t_->multiply(*t2.t_);  // t_ *= t2
		}

		void multiplyScalar(const ComplexOrRealType& scalar)
		{
			makeSureItsValid();
			t_->multiplyScalar(scalar);
		}

		const bool isEmpty() const { makeSureItsValid(); return t_->isEmpty(); }

		const bool metaEqual(const AlgebraType&) const
		{
			makeSureItsValid();
			return true;
		}

		static bool metaEqual(const AssignAndDestroy&, const AssignAndDestroy&)
		{
			return true;
		}

	private:

		AssignAndDestroy(const AssignAndDestroy&) = delete;

		AssignAndDestroy& operator=(const AssignAndDestroy&) = delete;

		void makeSureItsValid() const
		{
			if (isValid_) return;
			err("AssignAndDestroy: invalidated object!\n");
		}

		AlgebraType* t_;
		bool isValid_;
	};

	PairStringAuxType operator()(PsimagLite::String str, AuxiliaryType& aux) const
	{
		return PairStringAuxType(str, aux);
	}

	static bool isEmpty(const ResultType& term)
	{
		return term.isEmpty();
	}
};
}
#endif // SPECFORTARGETINGEXPRESSION_H
