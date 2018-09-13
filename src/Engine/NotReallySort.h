#ifndef NOT_REALLY_SORT_H
#define NOT_REALLY_SORT_H
#include "Qn.h"
#include "Vector.h"
#include <numeric>

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

		for (SizeType i = 0; i < n; ++i)
			hash[i] = computeHash(inQns[i], sizes, addOddToHash);

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

	static void computeSizes(VectorSizeType& sizes, const VectorQnType& inQns)
	{
		SizeType n = Qn::modalStruct.size(); // small number
		if (n == 0) return;
		sizes.resize(n);
		for (SizeType i = 0; i < n - 1; ++i)
			sizes[i] = findMaximum(inQns, i);
		sizes[n - 1] = 0; // should be unused
	}

	static SizeType findMaximum(const VectorQnType& inQns, SizeType index)
	{
		SizeType n = inQns.size(); // == very large number
		if (n == 0) return 1;
		SizeType max =  inQns[0].other[index];
		for (SizeType i = 1; i < n; ++i) {
			const SizeType val = inQns[i].other[index];
			if (val > max) max = val;
		}

		return max + 2;
	}
};
}
#endif // NOT_REALLY_SORT_H

