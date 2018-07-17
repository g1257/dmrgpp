#ifndef QN_H
#define QN_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "IndexOfItem.h"

namespace Dmrg {

class Qn {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<Qn>::Type VectorQnType;

	Qn(SizeType e, VectorSizeType szPlusConst, PairSizeType j, SizeType flavor)
	    : electrons(e), other(szPlusConst), jmPair(j), flavors(flavor)
	{}

	Qn(const Qn& q1, const Qn& q2)
	{
		electrons = q1.electrons + q2.electrons;
		SizeType n = q1.other.size();
		assert(q2.other.size() == n);
		other.resize(n);
		for (SizeType i = 0; i < n; ++i)
			other[i] = q1.other[i] + q2.other[i];

		jmPair.first = q1.jmPair.first + q2.jmPair.first;
		jmPair.second = q1.jmPair.second + q2.jmPair.second;
		flavors = q1.flavors; // ???
	}

	template<typename SomeInputType>
	void read(PsimagLite::String str, SomeInputType& io)
	{
		io.read(electrons, str + "/electrons");
		try {
			io.read(other, str + "/other");
		} catch (...) {}

		io.read(jmPair, str + "/jmPair");
		io.read(flavors, str + "/flavors");
	}

	void write(PsimagLite::String str, PsimagLite::IoNgSerializer& io) const
	{
		io.createGroup(str);
		io.write(str + "/electrons", electrons);
		io.write(str + "/other", other);
		io.write(str + "/jmPair", jmPair);
		io.write(str + "/flavors", flavors);
	}

	bool operator==(const Qn& a) const
	{
		return (a.electrons == electrons &&
		        vectorEqual(a.other) &&
		        pairEqual(a.jmPair) &&
		        flavors == a.flavors);
	}

	bool operator!=(const Qn& a) const
	{
		return !(*this == a);
	}

	void scale(SizeType sites,
	           SizeType totalSites,
	           ProgramGlobals::DirectionEnum direction,
	           bool isSu2)
	{
		Qn original = *this;
		SizeType mode = other.size();
		if (isSu2 && mode != 1)
			err("Qn::scale() expects mode==1 for SU(2)\n");

		if (direction == ProgramGlobals::INFINITE) {
			double flp = original.electrons*sites;
			electrons = static_cast<SizeType>(round(flp/totalSites));
			for (SizeType x = 0; x < mode; ++x) {
				flp = original.other[x]*sites;
				other[x] = static_cast<SizeType>(round(flp/totalSites));
			}

			flp = original.jmPair.first*sites;
			jmPair.first =  static_cast<SizeType>(round(flp/totalSites));
		}

		if (!isSu2) return;

		SizeType tmp =jmPair.first;
		PsimagLite::String str("SymmetryElectronsSz: FATAL: Impossible parameters ");
		bool flag = false;
		if (original.electrons & 1) {
			if (!(tmp & 1)) {
				flag = true;
				str += "electrons= " + ttos(original.electrons) + " is odd ";
				str += "and 2j= " +  ttos(tmp) + " is even.";
				tmp++;
			}
		} else {
			if (tmp & 1) {
				flag = true;
				str += "electrons= " + ttos(original.electrons) + " is even ";
				str += "and 2j= " +  ttos(tmp) + " is odd.";
				tmp++;
			}
		}

		if (flag && sites == totalSites) throw PsimagLite::RuntimeError(str);

		jmPair.first = tmp;
	}

	static void adjustQns(VectorQnType& outQns,
	                      const VectorSizeType& ints,
	                      SizeType mode)
	{
		SizeType n = ints.size();
		if (n == 0)
			err("adjustQns failed with n == 0\n");

		if (n % (mode + 1) != 0)
			err("adjustQns failed, n does not divide mode + 1\n");
		n /= (mode + 1);
		outQns.resize(n, Qn(0, VectorSizeType(mode), PairSizeType(0, 0), 0));
		for (SizeType i = 0; i < n; ++i) {
			assert(1 + i*(mode + 1) < ints.size());
			outQns[i].electrons = ints[1 + i*(mode + 1)];
			for (SizeType j = 0; j < mode; ++j) {
				SizeType k = (j == 0) ? 0 : j + 1;
				assert(k + i*mode < ints.size());
				assert(j < outQns[i].other.size());
				outQns[i].other[j] = ints[k + i*(mode + 1)];
			}
		}
	}

	static void notReallySort(VectorSizeType& outNumber,
	                          VectorQnType& outQns,
	                          VectorSizeType& partition,
	                          const VectorSizeType& inNumbers,
	                          const VectorQnType& inQns)
	{
		SizeType n = inNumbers.size();
		assert(n == inQns.size());

		outQns.clear();
		PsimagLite::Vector<VectorSizeType>::Type bucket;
		bucket.reserve(n);
		outQns.reserve(n);
		for (SizeType i = 0; i < n; ++i) {
			int x = PsimagLite::indexOfItemOrMinusOne(outQns, inQns[i]);
			if (x < 0) {
				outQns.push_back(inQns[i]);
				bucket.push_back(VectorSizeType(1, inNumbers[i]));
			} else {
				bucket[x].push_back(inNumbers[i]);
			}
		}

		SizeType buckets = bucket.size();
		SizeType counter = 0;
		outNumber.resize(n);
		partition.resize(buckets + 1);
		for (SizeType i = 0; i < buckets; ++i) {
			SizeType sizeOfThisBucket = bucket[i].size();
			partition[i] = counter;
			for (SizeType j = 0; j < sizeOfThisBucket; ++j) {
				assert(counter < outNumber.size());
				outNumber[counter++] = bucket[i][j];
			}
		}

		partition[buckets] = counter;
		assert(counter == outNumber.size());
	}

	static void qnToElectrons(VectorSizeType& electrons, const VectorQnType& qns)
	{
		electrons.resize(qns.size());
		for (SizeType i=0;i<qns.size();i++)
			electrons[i] = qns[i].electrons;
	}

	friend std::ostream& operator<<(std::ostream& os, const Qn& qn)
	{
		os<<"electrons="<<qn.electrons<<" ";
		os<<"other=";
		for (SizeType i = 0; i < qn.other.size(); ++i)
			os<<qn.other[i]<<",";
		os<<"  jmPair="<<qn.jmPair;
		return os;
	}

	SizeType electrons;
	VectorSizeType other;
	PairSizeType jmPair;
	SizeType flavors;

private:

	bool vectorEqual(const VectorSizeType& otherOther) const
	{
		SizeType n = otherOther.size();
		if (n != other.size()) return false;
		for (SizeType i = 0; i < n; ++i)
			if (otherOther[i] != other[i]) return false;
		return true;
	}

	bool pairEqual(const PairSizeType& otherJm) const
	{
		return (otherJm.first == jmPair.first &&
		        otherJm.second == jmPair.second);
	}
};

}
#endif // QN_H
