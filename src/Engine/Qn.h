#ifndef QN_H
#define QN_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Profiling.h"

namespace Dmrg {

class Qn {

public:

	enum ModalEnum { MODAL_SUM, MODAL_MODULO};

	struct ModalStruct {

		ModalStruct() : modalEnum(MODAL_SUM), extra(0) {}

		template<typename SomeInputType>
		void read(PsimagLite::String str, SomeInputType& io)
		{
			io.read(modalEnum, str + "/modalEnum");
			io.read(extra, str + "/extra");
		}

		void write(PsimagLite::String str, PsimagLite::IoNgSerializer& io) const
		{
			io.createGroup(str);
			io.write(str+ "/modalEnum", modalEnum);
			io.write(str + "/extra", extra);
		}

		ModalEnum modalEnum;
		short unsigned int extra;
	};

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<Qn>::Type VectorQnType;
	typedef PsimagLite::Vector<ModalStruct>::Type VectorModalStructType;

	Qn(bool odd, VectorSizeType szPlusConst, PairSizeType j, SizeType flavor)
	    : oddElectrons(odd), other(szPlusConst), jmPair(j), flavors(flavor)
	{
		if (modalStruct.size() == szPlusConst.size())
			return;

		if (szPlusConst.size() > 0)
			modalStruct.resize(szPlusConst.size());

		if (szPlusConst.size() == 0)
			err("Qn\n");
	}

	Qn(const Qn& q1, const Qn& q2)
	{
		oddElectrons = (q1.oddElectrons ^ q2.oddElectrons);
		SizeType n = q1.other.size();
		assert(q2.other.size() == n);
		assert(modalStruct.size() == n);

		other.resize(n);
		for (SizeType i = 0; i < n; ++i) {
			other[i] = q1.other[i] + q2.other[i];
			if (modalStruct[i].modalEnum == MODAL_MODULO)
				other[i] %= modalStruct[i].extra;
		}

		jmPair.first = q1.jmPair.first + q2.jmPair.first;
		jmPair.second = q1.jmPair.second + q2.jmPair.second;
		flavors = q1.flavors; // ???
	}

	template<typename SomeInputType>
	void read(PsimagLite::String str, SomeInputType& io)
	{
		io.read(oddElectrons, str + "/oddElectrons");
		try {
			io.read(other, str + "/other");
		} catch (...) {}

		io.read(jmPair, str + "/jmPair");
		io.read(flavors, str + "/flavors");

		if (modalStruct.size() == 0)
			io.read(modalStruct, "modalStruct");

		if (modalStruct.size() != other.size())
			err("Qn::read\n");
	}

	void write(PsimagLite::String str, PsimagLite::IoNgSerializer& io) const
	{
		static bool firstcall = true;
		if (firstcall) {
			io.write("modalStruct", modalStruct);
			firstcall = false;
		}

		io.createGroup(str);
		io.write(str + "/oddElectrons", oddElectrons);
		io.write(str + "/other", other);
		io.write(str + "/jmPair", jmPair);
		io.write(str + "/flavors", flavors);
	}

	bool operator==(const Qn& a) const
	{
		return (a.oddElectrons == oddElectrons &&
		        vectorEqualMaybeModal(a.other) &&
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
		if (isSu2 && mode != 2)
			err("Qn::scale() expects mode==1 for SU(2)\n");

		if (direction == ProgramGlobals::INFINITE) {
			double ts = totalSites;
			for (SizeType x = 0; x < mode; ++x) {
				double flp = original.other[x]*sites;
				other[x] = static_cast<SizeType>(round(flp/ts));
			}

			double flp = original.jmPair.first*sites;
			jmPair.first =  static_cast<SizeType>(round(flp/ts));
		}

		if (ifPresentOther0IsElectrons && other.size() > 0)
			oddElectrons = (other[0] & 1);

		if (!isSu2) return;

		assert(ifPresentOther0IsElectrons && other.size() > 0);

		SizeType tmp =jmPair.first;
		PsimagLite::String str("SymmetryElectronsSz: FATAL: Impossible parameters ");
		bool flag = false;
		if (original.oddElectrons) {
			if (!(tmp & 1)) {
				flag = true;
				str += "oddElectrons= " + ttos(original.oddElectrons) + " is odd ";
				str += "and 2j= " +  ttos(tmp) + " is even.";
				tmp++;
			}
		} else {
			if (tmp & 1) {
				flag = true;
				str += "oddElectrons= " + ttos(original.oddElectrons) + " is even ";
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
		SizeType modePlusOne = mode + 1;
		SizeType n = ints.size();
		if (n == 0)
			err("adjustQns failed with n == 0\n");

		if (n % modePlusOne != 0)
			err("adjustQns failed, n does not divide mode + 1\n");
		n /= modePlusOne;
		outQns.resize(n, Qn(false, VectorSizeType(modePlusOne), PairSizeType(0, 0), 0));

		for (SizeType i = 0; i < n; ++i) {
			assert(1 + i*modePlusOne < ints.size());
			SizeType tmp = ints[1 + i*modePlusOne];
			assert(outQns[i].other.size() > 0);
			outQns[i].other[0] = tmp;
			outQns[i].oddElectrons = (tmp & 1);
			for (SizeType j = 1; j < modePlusOne; ++j) {
				SizeType k = (j == 1) ? 0 : j;
				assert(k + i*modePlusOne < ints.size());
				assert(j < outQns[i].other.size());
				outQns[i].other[j] = ints[k + i*modePlusOne];
			}
		}
	}

	bool isDefinedOther() const
	{
		SizeType n = other.size();
		SizeType value = 1;
		SizeType total = sizeof(value)*8 - 1;
		value <<= total;
		for (SizeType i = 0; i < n; ++i)
			if (other[i] == value) return false;

		return true;
	}

	static void notReallySort(VectorSizeType& outNumber,
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

		outQns.clear();
		VectorSizeType count;
		count.reserve(n);
		VectorSizeType reverse(n);

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

	SizeType su2ElectronsBridge() const
	{
		assert(ifPresentOther0IsElectrons);
		assert(other.size() > 0);
		return other[0];
	}

	static void su2ElectronsBridge(VectorSizeType& v,
	                               const VectorQnType& qns)
	{
		SizeType n = qns.size();
		v.resize(n);
		for (SizeType i = 0; i < n; ++i)
			v[i] = qns[i].su2ElectronsBridge();
	}

	static Qn zero()
	{
		SizeType value = 1;
		SizeType total = sizeof(value)*8 - 1;
		value <<= total;
		return Qn(false, VectorSizeType(modalStruct.size(), value), PairSizeType(0, 0), 0);
	}

	friend std::ostream& operator<<(std::ostream& os, const Qn& qn)
	{
		os<<"oddElectrons="<<qn.oddElectrons<<" ";
		os<<"other=";
		for (SizeType i = 0; i < qn.other.size(); ++i)
			os<<qn.other[i]<<",";
		os<<"  jmPair="<<qn.jmPair;
		return os;
	}

	static VectorModalStructType modalStruct;
	static bool ifPresentOther0IsElectrons;
	bool oddElectrons;
	VectorSizeType other;
	PairSizeType jmPair;
	SizeType flavors;

private:

	bool vectorEqualMaybeModal(const VectorSizeType& otherOther) const
	{
		SizeType n = otherOther.size();
		if (n != other.size()) return false;
		assert(n == modalStruct.size());

		assert(isDefinedOther());

		for (SizeType i = 0; i < n; ++i) {
			if (modalStruct[i].modalEnum == MODAL_SUM) {
				if (otherOther[i] != other[i]) return false;
			} else {
				assert(modalStruct[i].modalEnum == MODAL_MODULO);
				assert(modalStruct[i].extra > 0);
				SizeType x = (otherOther[i] > other[i]) ? otherOther[i] - other[i] :
				                                          other[i] - otherOther[i];
				if ((x % modalStruct[i].extra) != 0) return false;
			}
		}

		return true;
	}

	bool pairEqual(const PairSizeType& otherJm) const
	{
		return (otherJm.first == jmPair.first &&
		        otherJm.second == jmPair.second);
	}

	// disable implicit conversion for 1st argument of ctor
	Qn(SizeType, VectorSizeType, PairSizeType, SizeType);
};
}
#endif // QN_H
