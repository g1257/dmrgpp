#ifndef TIMEEVOLVEFORTARGETINGEXPRESSION_H
#define TIMEEVOLVEFORTARGETINGEXPRESSION_H
#include "Vector.h"

namespace Dmrg {

template<typename PvectorsType>
class TimeEvolveForTargetingExpression {
	class OneTimeEvolution {

	public:

		typedef OneTimeEvolution ThisType;
		typedef typename PsimagLite::Vector<ThisType*>::Type VectorOneTimeEvolutionType;
		typedef typename PvectorsType::VectorWithOffsetType VectorWithOffsetType;
		typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
		typedef typename PvectorsType::RealType RealType;

		OneTimeEvolution(SizeType firstIndex,
		                 const VectorWithOffsetType& src,
		                 SizeType timeSteps,
		                 PvectorsType& pVectors)
		    : indices_(timeSteps), timesWithoutAdvancement_(0), time_(0)
		{
			indices_[0] = firstIndex;
			for (SizeType i = 1; i < timeSteps; ++i) {
				auto lambda = [this, i](SizeType ind) {
					indices_[i] = ind;
					return "|P" + ttos(ind) + ">";
				};

				pVectors.createNew(src, lambda);
			}
		}

		const VectorSizeType& indices() const { return indices_; }

		SizeType timesWithoutAdvancement() const
		{
			return timesWithoutAdvancement_;
		}

		RealType time() const { return time_; }

		void advanceTime(RealType tau)
		{
			time_ += tau;
		}

		void resetTimesWithoutAdvancement()
		{
			timesWithoutAdvancement_ = 1;
		}

		void incrementTimesWithoutAdvancement()
		{
			++timesWithoutAdvancement_;
		}

	private:

		VectorSizeType indices_;
		SizeType timesWithoutAdvancement_;
		RealType time_;
	};

public:

	typedef OneTimeEvolution OneTimeEvolutionType;
	typedef typename PsimagLite::Vector<OneTimeEvolution*>::Type VectorOneTimeEvolutionType;

	TimeEvolveForTargetingExpression()
	{}

	~TimeEvolveForTargetingExpression()
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i) {
			delete vEvolutions_[i];
			vEvolutions_[i] = nullptr;
		}
	}

	OneTimeEvolution* findThisEvolution(SizeType firstIndex) const
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i)
			if (vEvolutions_[i]->indices()[0] == firstIndex) return vEvolutions_[i];

		return nullptr;
	}

	void pushBack(OneTimeEvolution* ptr) { vEvolutions_.push_back(ptr); }

	SizeType size() const { return vEvolutions_.size(); }

private:

	VectorOneTimeEvolutionType vEvolutions_;

};
}
#endif // TIMEEVOLVEFORTARGETINGEXPRESSION_H
