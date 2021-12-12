#ifndef PVECTORS_H
#define PVECTORS_H
#include "InputNg.h"
#include "InputCheck.h"
#include "Pvector.h"
#include "TargetParamsTimeVectors.h"

namespace Dmrg {

template<typename TargetingBaseType>
class Pvectors {

public:

	typedef PsimagLite::InputNg<InputCheck>::Readable InputValidatorType;
	typedef typename TargetingBaseType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename TargetingBaseType::ApplyOperatorExpressionType ApplyOperatorExpressionType;
	typedef typename TargetingBaseType::ModelType ModelType;
	typedef typename ModelType::LeftRightSuperType LeftRightSuperType;
	typedef typename VectorWithOffsetType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef TargetParamsTimeVectors<ModelType> TargetParamsTimeVectorsType;
	typedef Pvector<typename VectorWithOffsetType::value_type> PvectorType;
	typedef typename PsimagLite::Vector<PvectorType*>::Type VectorPvectorType;
	typedef PsimagLite::Vector<bool>::Type VectorBoolType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	Pvectors(InputValidatorType& io,
	         const ApplyOperatorExpressionType& aoe,
	         const LeftRightSuperType& lrs)
	    : io_(io),
	      aoe_(aoe),
	      lrs_(lrs),
	      progress_("Pvectors"),
	      origPvectors_(0),
	      tstStruct_(io, "TargetingExpression", aoe_.model())
	{
		pvectorsFromInput(io);
	}

	~Pvectors()
	{
		for (SizeType i = 0; i < pVectors_.size(); ++i) {
			delete pVectors_[i];
			pVectors_[i] = 0;
		}
	}

	const PvectorType& operator()(SizeType ind) const
	{
		assert(ind < pVectors_.size());
		assert(pVectors_[ind]);
		return *pVectors_[ind];
	}

	ApplyOperatorExpressionType& aoeNonConst()
	{
		ApplyOperatorExpressionType* aoePtr = const_cast<ApplyOperatorExpressionType*>(&aoe_);
		return *aoePtr;
	}

	const ApplyOperatorExpressionType& aoe() const { return aoe_; }

	const LeftRightSuperType& lrs() const { return lrs_; }

	bool hasTimeEvolution(SizeType ind) const
	{
		assert(ind < pVectors_.size());
		assert(pVectors_[ind]);
		return pVectors_[ind]->hasTimeEvolution();
	}

	void setAsDone(SizeType ind)
	{
		assert(ind < pVectors_.size());
		assert(pVectors_[ind]);
		pVectors_[ind]->setAsDone();
	}

	void pushString(SizeType ind, PsimagLite::String str)
	{
		assert(ind < pVectors_.size());
		assert(pVectors_[ind]);
		pVectors_[ind]->pushString(str);
	}

	template<typename SomeLambdaType>
	void createNew(const VectorWithOffsetType& src, SomeLambdaType& lambda)
	{
		const SizeType ind = aoeNonConst().createPvector(src);
		const PsimagLite::String ename = lambda(ind);
		PvectorType* pnew = new PvectorType(ename);
		pnew->setAsDone();
		pVectors_.push_back(pnew);

		PsimagLite::OstringStream msgg(std::cout.precision());
		PsimagLite::OstringStream::OstringStreamType& msg = msgg();
		msg<<"P["<<ind<<"]="<<ename<<" created";
		progress_.printline(msgg, std::cout);
	}

	SizeType targets() const { return pVectors_.size(); }

	SizeType origPvectors() const { return origPvectors_; }

	int findInOrigNames(PsimagLite::String str) const
	{
		return findInNames(str, origPvectors_);
	}

	int findInAnyNames(PsimagLite::String str) const
	{
		return findInNames(str, pVectors_.size());
	}

	int findInNames(PsimagLite::String str, SizeType end) const
	{
		assert(end <= pVectors_.size());
		for (SizeType i = 0; i < end; ++i)
			if (pVectors_[i]->hasAnyName(str)) return i;

		return -1;
	}

	SizeType findPforThisExpression(PsimagLite::String str) const
	{
		const SizeType pvectors = pVectors_.size();
		for (SizeType i = 0; i < pvectors; ++i) {
			if (!pVectors_[i]->hasAnyName(str)) continue;
			return i;
		}

		throw PsimagLite::RuntimeError("findPforThisExpression: P not found\n");
	}

	void sumPvectors(SizeType ind0, SizeType ind1, PsimagLite::String p0PlusP1)
	{
		assert(ind0 < ind1);
		VectorWithOffsetType& v0 = aoeNonConst().targetVectorsNonConst(ind0);
		const VectorWithOffsetType& v1 =  aoe_.targetVectors(ind1);
		v0 += v1;
		pVectors_[ind0]->sum(*(pVectors_[ind1]), p0PlusP1);
	}

	void trimPvectors(PsimagLite::String str)
	{
		const SizeType tvs = aoe_.tvs();
		VectorBoolType used(tvs, false);
		assert(origPvectors_ <= tvs);
		for (SizeType i = 0; i < origPvectors_; ++i)
			used[i] = true;

		findUsedPvectors(used, str);

		for (SizeType i = 0; i < tvs; ++i) {
			if (used[i]) continue;
			aoeNonConst().destroyPvector(i);
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"P["<<i<<"] destroyed";
			progress_.printline(msgg, std::cout);
		}

		// this line is commented out because then the index of aoe.targetVectors
		// will be different than the index of pVectors
		// aoe.targetVectors is resized after all vectors have been computed
		// aoe().trimVectors();
		const SizeType tvsFinal = aoe_.tvs();
		if (tvs != tvsFinal) {
			PsimagLite::OstringStream msgg(std::cout.precision());
			PsimagLite::OstringStream::OstringStreamType& msg = msgg();
			msg<<"Number of target vectors is "<<tvsFinal<<" now";
			progress_.printline(msgg, std::cout);
		}
	}

	void initTimeVectors(SizeType timeSteps, RealType tau)
	{
		tstStruct_.times().resize(timeSteps);
		for (SizeType i = 0; i < timeSteps; ++i)
			tstStruct_.times()[i] = i*tau/(timeSteps - 1);

		ApplyOperatorExpressionType* aoePtr = const_cast<ApplyOperatorExpressionType*>(&aoe_);
		aoePtr->initTimeVectors(tstStruct_, io_);
	}

	const VectorWithOffsetType& getCurrentVectorConst(PsimagLite::String braOrKet) const
	{
		PsimagLite::GetBraOrKet getBraOrKet(braOrKet);
		if (getBraOrKet.isPvector()) {
			const SizeType pIndex = getBraOrKet.pIndex();
			if (pIndex >= aoe_.tvs())
				err("getVector: out of range for " + braOrKet + "\n");
			return aoe_.targetVectors(pIndex);
		} else if (getBraOrKet.isRvector()) {
			throw PsimagLite::RuntimeError("reserved vector\n");
		}

		const SizeType sectorIndex = getBraOrKet.sectorIndex();
		return *(aoe_.psiConst()[sectorIndex][getBraOrKet.levelIndex()]);
	}

private:

	void pvectorsFromInput(InputValidatorType& io)
	{
		checkObsolete(io);

		RealType gsWeight = 0;
		io.readline(gsWeight, "GsWeight=");

		const SizeType total = getPvectorsTotal(io);

		pVectors_.resize(total);
		PsimagLite::String tmp;
		RealType sum = 0.0;
		for (SizeType i = 0; i < total; ++i) {
			io.readline(tmp, "P" + ttos(i) + "=");
			PvectorType::checkSyntaxOfValue(tmp);
			pVectors_[i] = new PvectorType(tmp);
			sum += pVectors_[i]->weight();
		}

		if (sum == 0.0) return;

		RealType factor = (1.0 - gsWeight)/sum;
		for (SizeType i = 0; i < total; ++i)
			pVectors_[i]->multiplyWeight(factor);

		origPvectors_ = pVectors_.size();
	}

	static void checkObsolete(InputValidatorType& io)
	{
		bool hasObsolete = false;
		try {
			int tmp = 0;
			io.readline(tmp, "Pvectors=");
			hasObsolete = true;
		} catch (std::exception&) {}

		if (hasObsolete)
			err("Delete the Pvectors= line from the input; it's no longer needed\n");
	}

	static SizeType getPvectorsTotal(InputValidatorType& io)
	{
		VectorSizeType seen;
		for (auto it = io.map().begin(); it != io.map().end(); ++it) {
			int index = getPindex(it->first);
			if (index < 0) continue;
			seen.push_back(index);
		}

		const SizeType total = seen.size();
		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType iperm(total);
		sort.sort(seen, iperm);

		for (SizeType i = 0; i < total; ++i) {
			if (seen[i] == i) continue;
			err("Pvectors must not have holes\n");
		}

		return total;
	}

	static int getPindex(PsimagLite::String key)
	{
		const SizeType n = key.length();
		if (key.length() < 2) return -1;
		if (key[0] != 'P') return -1;

		PsimagLite::String buffer;
		for (SizeType i = 1; i < n; ++i) {
			if (key[i] < 48 || key[i] > 57) return -1;
			buffer += key[i];
		}

		return PsimagLite::atoi(buffer);
	}

	void findUsedPvectors(VectorBoolType& used, PsimagLite::String str) const
	{
		const SizeType len = str.size();
		for (SizeType i = 0; i < len; ++i) {

			if (str.substr(i, 2) != "|P") continue;

			PsimagLite::String buffer;
			SizeType j = i + 2;
			for (; j < len; ++j) {
				if (str[j] >= 48 && str[j] <= 57) {
					buffer += str[j];
				} else if (str[j] == '>') {
					break;
				} else {
					buffer = "";
					break;
				}
			}

			i = j;
			if (buffer == "") continue;
			const SizeType ind = PsimagLite::atoi(buffer);
			used[ind] = true;
		}
	}

	InputValidatorType& io_;
	const ApplyOperatorExpressionType& aoe_;
	const LeftRightSuperType& lrs_;
	PsimagLite::ProgressIndicator progress_;
	SizeType origPvectors_;
	VectorPvectorType pVectors_;
	TargetParamsTimeVectorsType tstStruct_;
};

}
#endif // PVECTORS_H
