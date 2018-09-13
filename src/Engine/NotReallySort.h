#ifndef NOT_REALLY_SORT_H
#define NOT_REALLY_SORT_H
#define USE_PTHREADS_OR_NOT_NG
#include "Qn.h"
#include "Vector.h"
#include <numeric>
#include "Sort.h"
#include "Concurrency.h"
#include "Parallelizer.h"

namespace Dmrg {

class NotReallySort {

public:

	typedef Qn::VectorQnType VectorQnType;
	typedef Qn::VectorSizeType VectorSizeType;

	// There exits objects Qn that have operator= (comparison)
	// notReallySort takes a vector of Qn objects as INPUT.1: inQns of size big
	// It also takes a vector of non-negative integers of size big as INPUT.2 : inNumbers
	// Let P be a permutation that "not really sorts" inQns, P is of size big
	// Let tmpQns be P applied to inQns; tmpQns[i] = inQns[P[i]], tmpQns is of size big
	// notReallySort fills the following vectors
	// OUTPUT.1: outNumber is P applied to inNumbers; outNumber[i] = inNumber[P[i]], of size big
	// OUTPUT.2: outQns[x] is the x-th unique tmpQns, of size small
	// OUTPUT.3: offset[x] = min {y; such that tmpQns[y] = outQns[x]}, of size small
	void operator()(VectorSizeType& outNumber,
	                VectorQnType& outQns,
	                VectorSizeType& offset,
	                const VectorSizeType& inNumbers,
	                const VectorQnType& inQns,
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
		firstPass(outQns, count, reverse, inQns, doNotSort);

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
			profiling->end("patches= " + ttos(numberOfPatches));
			delete profiling;
			profiling = 0;
		}
	}

private:

	class HashHelper {

	public:

		HashHelper(VectorSizeType& hash,
		           const VectorQnType& inQns,
		           const VectorSizeType& sizes,
		           bool addOdd)
		    : hash_(hash), inQns_(inQns), sizes_(sizes), addOdd_(addOdd)
		{}

		void doTask(SizeType taskNumber, SizeType)
		{
			hash_[taskNumber] = computeHash(inQns_[taskNumber], sizes_, addOdd_);
		}

		SizeType tasks() const { return inQns_.size(); }

	private:

		static SizeType computeHash(const Qn& qn, const VectorSizeType& sizes, bool addOdd)
		{
			SizeType n = qn.other.size(); // small number
			SizeType key = (addOdd && qn.oddElectrons) ? 1 : 0;
			SizeType scale = (addOdd) ? 2 : 1;

			for (SizeType i = 0; i < n; ++i) {
				key += qn.other[i]*scale;
				scale *= sizes[i];
			}

			return key;
		}

		VectorSizeType& hash_;
		const VectorQnType& inQns_;
		const VectorSizeType& sizes_;
		bool addOdd_;
	};

	static void firstPass(VectorQnType& outQns,
	                      VectorSizeType& count,
	                      VectorSizeType& reverse,
	                      const VectorQnType& inQns,
	                      bool doNotSort)
	{
		SizeType n = inQns.size();
		outQns.clear();
		if (n == 0) return;
		count.reserve(n);
		reverse.resize(n);

		VectorSizeType hash(n);
		VectorSizeType sizes;
		computeSizes(sizes, inQns);

		bool noNeedForOdd = (Qn::ifPresentOther0IsElectrons && Qn::modalStruct.size() > 0);
		bool hasAtLeastOneOdd = false;
		if (!noNeedForOdd) {
			for (SizeType i = 0; i < n; ++i) {
				hasAtLeastOneOdd = inQns[i].oddElectrons;
				if (hasAtLeastOneOdd) break;
			}
		}

		bool addOddToHash = (!noNeedForOdd && hasAtLeastOneOdd);

		SizeType threads = std::min(PsimagLite::Concurrency::codeSectionParams.npthreads, n);
		HashHelper helper(hash, inQns, sizes, addOddToHash);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<HashHelper> parallelizer(codeSectionParams);
		parallelizer.loopCreate(helper);

		PsimagLite::Sort<VectorSizeType> sort;
		VectorSizeType perm(n);

		if (doNotSort) {
			for (SizeType i = 0; i < n; ++i) perm[i] = i;
		} else {
			sort.sort(hash, perm);
		}

		SizeType j = 0;

		assert(n > 0);
		outQns.push_back(inQns[perm[0]]);
		count.push_back(1);
		reverse[perm[0]] = 0;
		for (SizeType i = 1; i < n; ++i) {
			const SizeType iperm = perm[i];
			if (hash[i - 1] == hash[i]) {
				++count[j];
				reverse[iperm] = j;
			} else {
				outQns.push_back(inQns[iperm]);
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

	class SizesHelper {

		typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

	public:

		SizesHelper(VectorSizeType& sizes, const VectorQnType& inQns, SizeType nthreads)
		    : sizes_(sizes), inQns_(inQns), otherSize_(Qn::modalStruct.size()), tmp_(nthreads)
		{
			for (SizeType thread = 0; thread < nthreads; ++thread)
				tmp_[thread].resize(otherSize_);
		}

		void doTask(SizeType taskNumber, SizeType thread)
		{
			for (SizeType index = 0; index < otherSize_ - 1; ++index) {

				const SizeType val = inQns_[taskNumber].other[index] + 2;

				// conditional assignment
				tmp_[thread][index] = (tmp_[thread][index] < val) ? val :
				                                                    tmp_[thread][index];
			}
		}

		SizeType tasks() const { return inQns_.size(); }

		void sync()
		{
			SizeType nthreads = tmp_.size();
			for (SizeType index = 0; index < otherSize_ - 1; ++index) {
				for (SizeType thread = 0; thread < nthreads; ++thread) {
					const SizeType val = tmp_[thread][index];
					sizes_[index] = (sizes_[index] < val) ? val : sizes_[index];
				}
			}

			sizes_[otherSize_ - 1] = 0; // should be unused
		}

	private:

		VectorSizeType& sizes_;
		const VectorQnType& inQns_;
		SizeType otherSize_;
		VectorVectorSizeType tmp_;
	};

	static void computeSizes(VectorSizeType& sizes, const VectorQnType& inQns)
	{
		const SizeType otherSize = Qn::modalStruct.size();
		if (otherSize == 0) return;
		SizeType n = inQns.size();
		sizes.resize(otherSize, 0);

		SizeType threads = std::min(PsimagLite::Concurrency::codeSectionParams.npthreads, n);
		SizesHelper helper(sizes, inQns, threads);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<SizesHelper> parallelizer(codeSectionParams);
		parallelizer.loopCreate(helper);
		helper.sync();
	}
};
}
#endif // NOT_REALLY_SORT_H

