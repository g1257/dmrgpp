#ifndef NOT_REALLY_SORT_H
#define NOT_REALLY_SORT_H
#include "Qn.h"
#include "Vector.h"
#include "Concurrency.h"

namespace Dmrg {

class NotReallySort {

	typedef Qn::VectorQnType VectorQnType;
	typedef Qn::VectorSizeType VectorSizeType;

	class NrsFirstPassHelper {

		typedef PsimagLite::Vector<VectorQnType>::Type VectorVectorQnType;
		typedef PsimagLite::Vector<VectorSizeType>::Type VectorVectorSizeType;

	public:

		NrsFirstPassHelper(const VectorQnType& inQns, SizeType threads)
		    : inQns_(inQns), count_(threads), outQns_(threads)
		{}

		void doTask(SizeType taskNumber, SizeType threadNum)
		{
			int x = PsimagLite::indexOrMinusOne(outQns_[threadNum], inQns_[taskNumber]);
			if (x < 0) {
				outQns_[threadNum].push_back(inQns_[taskNumber]);
				count_[threadNum].push_back(1);
			} else {
				++count_[threadNum][x];
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
					count[x] = count_[t][i];
					inQnsReduced[x++] = outQns_[t][i];
				}

				count_[t].clear();
				outQns_[t].clear();
			}

			count_.clear();
			outQns_.clear();

			outQns.clear();
			count.reserve(nreduced);
			SizeType n = inQns_.size();
			reverse.resize(n);

			SizeType bigI = 0;
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

				for (SizeType j = 0; j < c; ++j)
					reverse[bigI++] = x;
			}
		}

	private:

		const VectorQnType& inQns_;
		VectorVectorSizeType count_;
		VectorVectorQnType outQns_;
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
		nrsFirstPass(outQns, count, reverse, inQns);

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

	static void nrsFirstPass(VectorQnType& outQns,
	                         VectorSizeType& count,
	                         VectorSizeType& reverse,
	                         const VectorQnType& inQns)
	{
#ifndef USE_PTHREADS
		return nrsFirstPassSerial(outQns, count, reverse, inQns);
#endif

		SizeType threads = PsimagLite::Concurrency::codeSectionParams.npthreads;
		NrsFirstPassHelper helper(inQns, threads);
		PsimagLite::CodeSectionParams codeSectionParams(threads);
		PsimagLite::Parallelizer<NrsFirstPassHelper> parallelizer(codeSectionParams);
		parallelizer.loopCreate(helper);
		helper.sync(outQns, count, reverse);
	}

	static void nrsFirstPassSerial(VectorQnType& outQns,
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
