#ifndef NOT_REALLY_SORT_H
#define NOT_REALLY_SORT_H
#include "Qn.h"
#include "Vector.h"
#include "Concurrency.h"

namespace Dmrg {

class NotReallySort {

	typedef Qn::VectorQnType VectorQnType;
	typedef Qn::VectorSizeType VectorSizeType;

	class FirstPassHelper {

		typedef PsimagLite::Vector<VectorQnType>::Type VectorVectorQnType;
		typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

	public:

		FirstPassHelper(const VectorQnType& inQns, SizeType threads)
		    : inQns_(inQns), count_(threads), outQns_(threads), reverse_(inQns.size())
		{}

		void doTask(SizeType taskNumber, SizeType threadNum)
		{
			int x = PsimagLite::indexOrMinusOne(outQns_[threadNum], inQns_[taskNumber]);
			if (x < 0) {
				outQns_[threadNum].push_back(inQns_[taskNumber]);
				count_[threadNum].push_back(1);
				reverse_[taskNumber] = count_[threadNum].size() - 1;
			} else {
				++count_[threadNum][x];
				reverse_[taskNumber] = x;
			}
		}

		SizeType tasks() const { return inQns_.size(); }

		void sync(VectorQnType& outQns,
		          VectorSizeType& count,
		          VectorSizeType& reverse)
		{
			SizeType threads = count_.size();
			SizeType nreduced = 0;
			for (SizeType t = 0; t < threads; ++t)
				nreduced += count_[t].size();

			VectorSizeType countReduced(nreduced);
			VectorQnType inQnsReduced(nreduced, Qn::zero());

			SizeType x = 0;
			for (SizeType t = 0; t < threads; ++t) {
				SizeType m = count_[t].size();
				for (SizeType i = 0; i < m; ++i) {
					countReduced[x] = count_[t][i];
					inQnsReduced[x++] = outQns_[t][i];
				}
			}

			count.reserve(nreduced);

			for (SizeType ii = 0; ii < nreduced; ++ii) {
				int x = PsimagLite::indexOrMinusOne(outQns, inQnsReduced[ii]);
				SizeType c = countReduced[ii];
				if (x < 0) {
					outQns.push_back(inQnsReduced[ii]);
					count.push_back(c);
					x = count.size() - 1;
				} else {
					count[x] += c;
				}
			}

			reverse = reverse_; // <--- THIS IS WRONG!
		}

	private:

		const VectorQnType& inQns_;
		VectorVectorSizeType count_;
		VectorVectorQnType outQns_;
		VectorSizeType reverse_;
	};

public:

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
		firstPass(outQns, count, reverse, inQns);

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

	static void firstPass(VectorQnType& outQns,
	                      VectorSizeType& count,
	                      VectorSizeType& reverse,
	                      const VectorQnType& inQns)
	{
		return firstPassSerial(outQns, count, reverse, inQns);

		/* SizeType threads = PsimagLite::Concurrency::codeSectionParams.npthreads;
		FirstPassHelper helper(inQns, threads);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<FirstPassHelper> parallelizer(codeSectionParams);
		parallelizer.loopCreate(helper);
		helper.sync(outQns, count, reverse);*/
	}

	/*static void firstPass(VectorQnType& outQns,
	                      VectorSizeType& count,
	                      VectorSizeType& reverse,
	                      const VectorQnType& inQns)
	{
		VectorQnType outQns2;
		VectorSizeType  count2;
		VectorSizeType reverse2;
		firstPassSerial(outQns2, count2, reverse2, inQns);

		SizeType threads = PsimagLite::Concurrency::codeSectionParams.npthreads;
		FirstPassHelper helper(inQns, threads);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<FirstPassHelper> parallelizer(codeSectionParams);
		parallelizer.loopCreate(helper);
		helper.sync(outQns, count, reverse);
		assert(compare(outQns, outQns2));
		assert(compare(count, count2));
		assert(compare(reverse, reverse2));
	}

	template<typename T>
	static bool compare(const T& v1, const T& v2)
	{
		if (v1.size() != v2.size()) return false;
		for (SizeType i = 0; i < v1.size(); ++i)
			if (v1[i] != v2[i]) return false;
		return true;
	}*/

	static void firstPassSerial(VectorQnType& outQns,
	                            VectorSizeType& count,
	                            VectorSizeType& reverse,
	                            const VectorQnType& inQns)
	{
		SizeType n = inQns.size();
		outQns.clear();
		count.reserve(n);
		reverse.resize(n);

		// 1^st pass over data
		for (SizeType i = 0; i < n; ++i) {
			int x = PsimagLite::indexOrMinusOne(outQns, inQns[i]);
			if (x < 0) {
				outQns.push_back(inQns[i]);
				count.push_back(1);
				reverse[i] = count.size() - 1;
			} else {
				++count[x];
				reverse[i] = x;
			}
		}
	}
};
}
#endif // NOT_REALLY_SORT_H
