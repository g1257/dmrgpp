#ifndef NOT_REALLY_SORT_H
#define NOT_REALLY_SORT_H
#define USE_PTHREADS_OR_NOT_NG
#include "Qn.h"
#include "Vector.h"
#include <numeric>
#include "Sort.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include <tr1/unordered_map>
#include "PairOfQns.h"
#include "Array.h"

namespace std {

namespace tr1 {

template<>
class hash<Dmrg::Qn> {

	typedef Dmrg::Qn::VectorQnType VectorQnType;
	typedef Dmrg::Array<Dmrg::PairOfQns> VectorLikeQnType;
	typedef Dmrg::Qn::VectorSizeType VectorSizeType;

public:

	hash(VectorSizeType& hash,
	     const VectorQnType& inQns,
	     bool addOdd)
	    : hash_(hash), inQns_(&inQns), inQns2_(0), addOdd_(addOdd)
	{}

	hash(VectorSizeType& hash,
	     const VectorLikeQnType& inQns,
	     bool addOdd)
	    : hash_(hash), inQns_(0), inQns2_(&inQns), addOdd_(addOdd)
	{}

	void doTask(SizeType taskNumber, SizeType)
	{
		hash_[taskNumber] = (inQns2_) ? operator()((*inQns2_)[taskNumber]) :
		                                operator()((*inQns_)[taskNumber]);
	}

	SizeType tasks() const
	{
		return (inQns2_) ? inQns2_->size() : inQns_->size();
	}

	SizeType operator()(const Dmrg::PairOfQns& qnPair) const
	{
		return operator()(qnPair.make());
	}

	SizeType operator()(const Dmrg::Qn& qn) const
	{
		const SizeType offset = 8; // 8 bits
		const SizeType n = qn.other.size(); // small number
		SizeType key = (addOdd_ && qn.oddElectrons) ? 1 : 0;
		SizeType bits = (addOdd_) ? 1 : 0;

		for (SizeType i = 0; i < n; ++i) {
			SizeType val = qn.other[i];
			val <<= bits;
			key += val;
			bits += offset;
		}

		return key;
	}

private:

	VectorSizeType& hash_;
	const VectorQnType* inQns_;
	const VectorLikeQnType* inQns2_;
	bool addOdd_;
};
}
}

namespace Dmrg {

class NotReallySort {

public:

	typedef Qn::VectorQnType VectorQnType;
	typedef Qn::VectorSizeType VectorSizeType;

	enum AlgoEnum {ALGO_UMAP, ALGO_CUSTOM};

	NotReallySort()
	{
		if (ProgramGlobals::notReallySortAlgo == "custom")
			algo_ = ALGO_CUSTOM;
		else
			algo_ = ALGO_UMAP;
	}

	// There exits objects Qn that have operator= (comparison)
	// notReallySort takes a vector of Qn objects as INPUT.1: inQns of size big
	// It also takes a vector of non-negative integers of size big as INPUT.2 : inNumbers
	// Let P be a permutation that "not really sorts" inQns, P is of size big
	// Let tmpQns be P applied to inQns; tmpQns[i] = inQns[P[i]], tmpQns is of size big
	// notReallySort fills the following vectors
	// OUTPUT.1: outNumber is P applied to inNumbers; outNumber[i] = inNumber[P[i]], of size big
	// OUTPUT.2: outQns[x] is the x-th unique tmpQns, of size small
	// OUTPUT.3: offset[x] = min {y; such that tmpQns[y] = outQns[x]}, of size small
	template<typename SomeVectorLikeQnType>
	void operator()(VectorSizeType& outNumber,
	                VectorQnType& outQns,
	                VectorSizeType& offset,
	                const VectorSizeType& inNumbers,
	                const SomeVectorLikeQnType& inQns,
	                bool doNotSort,
	                ProgramGlobals::VerboseEnum verbose)
	{
		SizeType n = inNumbers.size();
		assert(n == inQns.size());
		PsimagLite::Profiling* profiling = (verbose) ? new PsimagLite::Profiling("notReallySort",
		                                                                         "n= " + ttos(n),
		                                                                         std::cout) : 0;

		VectorSizeType count;
		VectorSizeType reverse;

		// 1^st pass over data
		if (algo_ == ALGO_CUSTOM)
			firstPass(outQns, count, reverse, inQns, doNotSort);
		else
			firstPass2(outQns, count, reverse, inQns, doNotSort);

		// perform prefix sum
		SizeType numberOfPatches = count.size();
		offset.resize(numberOfPatches + 1);
		offset[0] = 0;
		for (SizeType ipatch = 0; ipatch < numberOfPatches; ++ipatch)
			offset[ipatch + 1] = offset[ipatch] + count[ipatch];

		// 2^nd pass over data
		outNumber.resize(n);
		std::fill(count.begin(), count.end(), 0);
		for (SizeType i = 0; i < n; ++i) {
			SizeType x = reverse[i];
			assert(x < offset.size() && x < count.size());
			SizeType outIndex = offset[x] + count[x];
			outNumber[outIndex] = inNumbers[i];
			++count[x];
		}

		if (profiling) {
			profiling->end("patches= " + ttos(numberOfPatches) +
			               " bitwise key, algo=" + ProgramGlobals::notReallySortAlgo);
			delete profiling;
			profiling = 0;
		}
	}

private:

	template<typename SomeVectorLikeQnType>
	void firstPass(VectorQnType& outQns,
	               VectorSizeType& count,
	               VectorSizeType& reverse,
	               const SomeVectorLikeQnType& inQns,
	               bool doNotSort)
	{
		SizeType n = inQns.size();
		outQns.clear();
		if (n == 0) return;
		count.reserve(n);
		reverse.resize(n);

		VectorSizeType hash(n);

		bool noNeedForOdd = (Qn::ifPresentOther0IsElectrons && Qn::modalStruct.size() > 0);
		bool hasAtLeastOneOdd = false;
		if (!noNeedForOdd) {
			for (SizeType i = 0; i < n; ++i) {
				hasAtLeastOneOdd = makeQnIfNeeded(inQns[i]).oddElectrons;
				if (hasAtLeastOneOdd) break;
			}
		}

		bool addOddToHash = (!noNeedForOdd && hasAtLeastOneOdd);

		SizeType threads = std::min(PsimagLite::Concurrency::codeSectionParams.npthreads, n);
		std::tr1::hash<Qn> helper(hash, inQns, addOddToHash);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<std::tr1::hash<Qn> > parallelizer(codeSectionParams);
		parallelizer.loopCreate(helper);

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType perm(n);

		if (doNotSort) {
			checkThatHashIsNotReallySorted(hash);
			for (SizeType i = 0; i < n; ++i) perm[i] = i;
		} else {
			sort.sort(hash, perm);
		}

		SizeType j = 0;

		assert(n > 0);
		outQns.push_back(makeQnIfNeeded(inQns[perm[0]]));
		count.push_back(1);
		reverse[perm[0]] = 0;
		for (SizeType i = 1; i < n; ++i) {
			const SizeType iperm = perm[i];
			if (hash[i - 1] == hash[i]) {
				++count[j];
				reverse[iperm] = j;
			} else {
				outQns.push_back(makeQnIfNeeded(inQns[iperm]));
				count.push_back(1);
				++j;
				// assert(j == count.size() - 1);
				// assert(count.size() == outQns.size());
				reverse[iperm] = j;
			}
		}

		//checkSum(count, n);
		//checkReverse(inQns, reverse, outQns);
	}

	// only for debugging
	static void checkThatHashIsNotReallySorted(const VectorSizeType& h)
	{
#ifdef NDEBUG
		return;
#else

		SizeType n = h.size();
		if (n == 0) return;

		VectorSizeType seen(1, h[0]);
		for (SizeType i = 1; i < n; ++i) {
			if (h[i] == h[i - 1])
				continue;

			assert(PsimagLite::indexOrMinusOne(seen, h[i]) < 0);
			seen.push_back(h[i]);
		}
#endif
	}

	template<typename SomeVectorLikeQnType>
	void firstPass2(VectorQnType& outQns,
	                VectorSizeType& count,
	                VectorSizeType& reverse,
	                const SomeVectorLikeQnType& inQns,
	                bool doNotSort)
	{
		if (doNotSort)
			return firstPass(outQns, count, reverse, inQns, doNotSort);

		SizeType n = inQns.size();
		outQns.clear();
		if (n == 0) return;
		count.reserve(n);
		reverse.reserve(n);

		VectorSizeType hash(1);

		bool noNeedForOdd = (Qn::ifPresentOther0IsElectrons && Qn::modalStruct.size() > 0);
		bool hasAtLeastOneOdd = false;
		if (!noNeedForOdd) {
			for (SizeType i = 0; i < n; ++i) {
				hasAtLeastOneOdd = makeQnIfNeeded(inQns[i]).oddElectrons;
				if (hasAtLeastOneOdd) break;
			}
		}

		bool addOddToHash = (!noNeedForOdd && hasAtLeastOneOdd);

		std::tr1::hash<Qn> helper(hash, inQns, addOddToHash);
		std::tr1::unordered_map<Qn, SizeType> umap(10, helper);

		for (SizeType i = 0; i < n; ++i) {
			const Qn qn = makeQnIfNeeded((inQns[i]));
			const  bool isSeenBefore =  (umap.count(qn) > 0);

			if (!isSeenBefore) {
				outQns.push_back(qn);
				count.push_back(1);
				reverse.push_back(count.size() - 1);
				umap[qn] = reverse[i];
			} else {
				const SizeType x = umap[qn];
				++count[x];
				reverse.push_back(x);
			}
		}
	}

	static Qn makeQnIfNeeded(const Qn& qn)
	{
		return qn;
	}

	static Qn makeQnIfNeeded(const PairOfQns& qnPair)
	{
		return qnPair.make();
	}

	AlgoEnum algo_;
};
}
#endif // NOT_REALLY_SORT_H

