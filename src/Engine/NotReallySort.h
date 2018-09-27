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

public:

	typedef Dmrg::Qn::VectorQnType VectorQnType;
	typedef Dmrg::Qn::VectorSizeType VectorSizeType;

	hash(VectorSizeType& hash,
	     const VectorQnType& inQns,
	     bool addOdd)
	    : hash_(hash), inQns_(inQns), addOdd_(addOdd)
	{}

	void doTask(SizeType taskNumber, SizeType)
	{
		assert(inQns_.size() > taskNumber);
		hash_[taskNumber] = operator()(inQns_[taskNumber]);
	}

	SizeType tasks() const
	{
		return inQns_.size();
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
	const VectorQnType& inQns_;
	bool addOdd_;
};

template<>
class hash<Dmrg::PairOfQns> {

public:

	typedef Dmrg::Array<Dmrg::PairOfQns> VectorLikeQnType;
	typedef Dmrg::Qn::VectorSizeType VectorSizeType;


	hash(VectorSizeType& hash,
	     const VectorLikeQnType& inQns,
	     bool addOdd)
	    : hash_(hash), inQns_(inQns), addOdd_(addOdd)
	{}

	SizeType operator()(const Dmrg::PairOfQns& qnPair) const
	{
		return qnPair.hash(addOdd_);
	}

	void doTask(SizeType taskNumber, SizeType)
	{
		assert(inQns_.size() > taskNumber);
		hash_[taskNumber] = operator()(inQns_[taskNumber]);
	}

	SizeType tasks() const
	{
		return inQns_.size();
	}

private:

	VectorSizeType& hash_;
	const VectorLikeQnType& inQns_;
	bool addOdd_;
};
}
}

namespace Dmrg {

class NotReallySort {

public:

	typedef Qn::VectorQnType VectorQnType;
	typedef Qn::VectorSizeType VectorSizeType;
	typedef typename std::tr1::hash<Dmrg::PairOfQns>::VectorLikeQnType VectorLikeQnType;

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
	                SizeType initialSizeOfHashTable,
	                ProgramGlobals::VerboseEnum verbose)
	{
		typedef typename SomeVectorLikeQnType::value_type PairOfQnsOrJustQnType;

		SizeType n = inNumbers.size();
		assert(n == inQns.size());
		PsimagLite::Profiling* profiling = (verbose) ? new PsimagLite::Profiling("notReallySort",
		                                                                         "n= " + ttos(n),
		                                                                         std::cout) : 0;

		VectorSizeType count;

		if (algo_ == ALGO_CUSTOM || doNotSort) {
			VectorSizeType reverse;
			firstPassCustom(outQns, count, reverse, inQns, doNotSort);
			secondPassCustom(outNumber, offset, count, reverse, inNumbers, inQns);
		} else {
			VectorSizeType hash(1);

			const bool noNeedForOdd = (Qn::ifPresentOther0IsElectrons &&
			                           Qn::modalStruct.size() > 0);
			const bool hasAtLeastOneOdd = (noNeedForOdd) ? false : getOddElectrons(inQns);
			const bool addOddToHash = (!noNeedForOdd && hasAtLeastOneOdd);

			std::tr1::hash<PairOfQnsOrJustQnType> helper(hash, inQns, addOddToHash);
			std::tr1::unordered_map<PairOfQnsOrJustQnType, SizeType> umap(initialSizeOfHashTable,
			                                                              helper);
			firstPassUmap(outQns, count, umap, inQns);
			secondPassUmap(outNumber, offset, count, umap, inNumbers, inQns);
		}

		SizeType numberOfPatches = count.size();

		if (profiling) {
			profiling->end("patches= " + ttos(numberOfPatches) +
			               " bitwise key, algo=" + ProgramGlobals::notReallySortAlgo);
			delete profiling;
			profiling = 0;
		}
	}

private:

	template<typename SomeVectorLikeQnType>
	void firstPassCustom(VectorQnType& outQns,
	                     VectorSizeType& count,
	                     VectorSizeType& reverse,
	                     const SomeVectorLikeQnType& inQns,
	                     bool doNotSort)
	{
		typedef typename SomeVectorLikeQnType::value_type PairOfQnsOrJustQnType;

		SizeType n = inQns.size();
		outQns.clear();
		if (n == 0) return;
		count.reserve(n);
		reverse.resize(n);

		VectorSizeType hash(n);

		const bool noNeedForOdd = (Qn::ifPresentOther0IsElectrons &&
		                           Qn::modalStruct.size() > 0);
		const bool hasAtLeastOneOdd = (noNeedForOdd) ? false : getOddElectrons(inQns);
		const bool addOddToHash = (!noNeedForOdd && hasAtLeastOneOdd);

		SizeType threads = std::min(PsimagLite::Concurrency::codeSectionParams.npthreads, n);
		std::tr1::hash<PairOfQnsOrJustQnType> helper(hash, inQns, addOddToHash);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<std::tr1::hash<PairOfQnsOrJustQnType> >
		        parallelizer(codeSectionParams);
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

	template<typename SomeVectorLikeQnType>
	void secondPassCustom(VectorSizeType& outNumber,
	                      VectorSizeType& offset,
	                      VectorSizeType& count,
	                      const VectorSizeType& reverse,
	                      const VectorSizeType& inNumbers,
	                      const SomeVectorLikeQnType& inQns)
	{
		SizeType n = inNumbers.size();
		assert(n == inQns.size());

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
	void firstPassUmap(VectorQnType& outQns,
	                   VectorSizeType& count,
	                   std::tr1::unordered_map<typename SomeVectorLikeQnType::value_type,
	                   SizeType>& umap,
	                   const SomeVectorLikeQnType& inQns)
	{
		typedef typename SomeVectorLikeQnType::value_type PairOfQnsOrJustQnType;

		SizeType n = inQns.size();
		outQns.clear();
		if (n == 0) return;
		count.reserve(n);

		for (SizeType i = 0; i < n; ++i) {
			const PairOfQnsOrJustQnType& qn = inQns[i];
			const bool isSeenBefore =  (umap.count(qn) > 0);

			if (!isSeenBefore) {
				outQns.push_back(makeQnIfNeeded(qn));
				count.push_back(1);
				umap[qn] = count.size() - 1;
			} else {
				const SizeType x = umap[qn];
				++count[x];
			}
		}
	}

	template<typename SomeVectorLikeQnType>
	void secondPassUmap(VectorSizeType& outNumber,
	                    VectorSizeType& offset,
	                    VectorSizeType& count,
	                    std::tr1::unordered_map<typename SomeVectorLikeQnType::value_type,
	                    SizeType>& umap,
	                    const VectorSizeType& inNumbers,
	                    const SomeVectorLikeQnType& inQns)
	{
		SizeType n = inNumbers.size();
		assert(n == inQns.size());

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
			SizeType x = umap[inQns[i]];
			assert(x < offset.size() && x < count.size());
			SizeType outIndex = offset[x] + count[x];
			outNumber[outIndex] = inNumbers[i];
			++count[x];
		}
	}

	static bool getOddElectrons(const VectorQnType& inQns)
	{
		bool hasAtLeastOneOdd = false;
		SizeType n = inQns.size();
		for (SizeType i = 0; i < n; ++i) {
			hasAtLeastOneOdd = inQns[i].oddElectrons;
			if (hasAtLeastOneOdd) break;
		}

		return hasAtLeastOneOdd;
	}

	static bool getOddElectrons(const Array<PairOfQns>& inQns)
	{
		bool hasAtLeastOneOdd = false;
		SizeType n = inQns.size();
		for (SizeType i = 0; i < n; ++i) {
			hasAtLeastOneOdd = inQns[i].oddElectrons();
			if (hasAtLeastOneOdd) break;
		}

		return hasAtLeastOneOdd;
	}

	static const Qn& makeQnIfNeeded(const Qn& qn) { return qn; }

	static Qn makeQnIfNeeded(const PairOfQns& qn) { return qn.make(); }

	AlgoEnum algo_;
};
}
#endif // NOT_REALLY_SORT_H

