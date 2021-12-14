#ifndef GROUP_OF_ONE_TIME_EVOLUTIONS_H
#define GROUP_OF_ONE_TIME_EVOLUTIONS_H
#include "Vector.h"
#include "GetBraOrKet.h"

namespace Dmrg {

template<typename PvectorsType>
class GroupOfOneTimeEvolutions {
	class OneTimeEvolution {

	public:

		typedef OneTimeEvolution ThisType;
		typedef typename PsimagLite::Vector<ThisType*>::Type VectorOneTimeEvolutionType;
		typedef typename PvectorsType::VectorWithOffsetType VectorWithOffsetType;
		typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
		typedef typename PvectorsType::RealType RealType;

		OneTimeEvolution(SizeType firstIndex,
		                 const VectorWithOffsetType& src,
		                 PsimagLite::String srcKet,
		                 SizeType disposition,
		                 SizeType timeSteps,
		                 PvectorsType& pVectors)
		    : indices_(timeSteps),
		      srcKet_(srcKet),
		      disposition_(disposition),
		      timesWithoutAdvancement_(0),
		      time_(0)
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

		bool hasKet(SizeType ket) const
		{
			const SizeType n = indices_.size();
			for (SizeType i = 0; i < n; ++i)
				if (indices_[i] == ket) return true;

			return false;
		}

		int sourceToDestroy() const
		{
			return (disposition_ != 2) ? -1 : getPindexOrMinusOne(srcKet_);
		}

	private:

		static int getPindexOrMinusOne(PsimagLite::String braOrKet)
		{
			PsimagLite::GetBraOrKet getBraOrKet(braOrKet);
			return (getBraOrKet.isPvector()) ? getBraOrKet.pIndex() : -1;
		}

		VectorSizeType indices_;
		PsimagLite::String srcKet_;
		SizeType disposition_;
		SizeType timesWithoutAdvancement_;
		RealType time_;
	};

public:

	typedef OneTimeEvolution OneTimeEvolutionType;
	typedef typename PsimagLite::Vector<OneTimeEvolution*>::Type VectorOneTimeEvolutionType;
	typedef typename PvectorsType::RealType RealType;

	GroupOfOneTimeEvolutions() {}

	GroupOfOneTimeEvolutions(const GroupOfOneTimeEvolutions&) = delete;

	GroupOfOneTimeEvolutions& operator=(const GroupOfOneTimeEvolutions&) = delete;

	~GroupOfOneTimeEvolutions()
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

	RealType getTimeForKet(SizeType ket) const
	{
		const SizeType n = vEvolutions_.size();
		for (SizeType i = 0; i < n; ++i) {
			if (vEvolutions_[i]->hasKet(ket)) return vEvolutions_[i]->time();
		}

		return 0;
	}

private:

	VectorOneTimeEvolutionType vEvolutions_;

};
}
#endif // GROUP_OF_ONE_TIME_EVOLUTIONS_H
