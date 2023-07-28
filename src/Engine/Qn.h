#ifndef QN_H
#define QN_H
#include "Array.h"
#include "Io/IoNg.h"
#include "Profiling.h"
#include "ProgramGlobals.h"
#include "Vector.h"

namespace Dmrg
{

class Qn
{

public:

	enum ModalEnum { MODAL_SUM,
		MODAL_MODULO };

	struct ModalStruct {

		ModalStruct()
		    : modalEnum(MODAL_SUM)
		    , extra(0)
		{
		}

		template <typename SomeInputType>
		void read(PsimagLite::String str, SomeInputType& io)
		{
			io.read(modalEnum, str + "/modalEnum");
			io.read(extra, str + "/extra");
		}

		void write(PsimagLite::String str,
		    PsimagLite::IoNgSerializer& io,
		    typename PsimagLite::IoNgSerializer::WriteMode wM = PsimagLite::IoNgSerializer::NO_OVERWRITE) const
		{
			if (wM != PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
				io.createGroup(str);
			io.write(str + "/modalEnum", modalEnum, wM);
			io.write(str + "/extra", extra, wM);
		}

		ModalEnum modalEnum;
		short unsigned int extra;
	};

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef std::pair<SizeType, SizeType> PairSizeType;
	typedef PsimagLite::Vector<Qn>::Type VectorQnType;
	typedef PsimagLite::Vector<ModalStruct>::Type VectorModalStructType;

	Qn(bool odd, VectorSizeType szPlusConst, PairSizeType j, SizeType flavor)
	    : oddElectrons(odd)
	    , other(szPlusConst)
	    , jmPair(j)
	    , flavors(flavor)
	{
		if (modalStruct.size() == szPlusConst.size()) {
			modularize();
			return;
		}

		if (szPlusConst.size() > 0) {
			modalStruct.resize(szPlusConst.size());
			modularize();
		}

		if (szPlusConst.size() == 0)
			modalStruct.clear();
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

	template <typename SomeInputType>
	void read(PsimagLite::String str, SomeInputType& io)
	{
		io.read(oddElectrons, str + "/oddElectrons");
		try {
			VectorSizeType otherVector;
			io.read(otherVector, str + "/other");
			other.fromStdVector(otherVector);
		} catch (...) {
		}

		io.read(jmPair, str + "/jmPair");
		io.read(flavors, str + "/flavors");

		if (modalStruct.size() == 0)
			io.read(modalStruct, "modalStruct");

		if (modalStruct.size() != other.size())
			err("Qn::read\n");
	}

	void write(PsimagLite::String str,
	    PsimagLite::IoNgSerializer& io,
	    typename PsimagLite::IoNgSerializer::WriteMode wM = PsimagLite::IoNgSerializer::NO_OVERWRITE) const
	{
		try {
			io.read(modalStruct, "modalStruct");
		} catch (...) {
			io.write("modalStruct", modalStruct);
		}

		if (wM != PsimagLite::IoNgSerializer::ALLOW_OVERWRITE)
			io.createGroup(str);
		io.write(str + "/oddElectrons", oddElectrons, wM);
		VectorSizeType otherVector;
		other.toStdVector(otherVector);
		io.write(str + "/other", otherVector, wM);
		io.write(str + "/jmPair", jmPair, wM);
		io.write(str + "/flavors", flavors, wM);
	}

	void overwrite(PsimagLite::String str, PsimagLite::IoNgSerializer& io) const
	{
		const PsimagLite::IoNgSerializer::WriteMode mode = PsimagLite::IoNgSerializer::ALLOW_OVERWRITE;
		io.write(str + "/oddElectrons", oddElectrons, mode);
		VectorSizeType otherVector;
		other.toStdVector(otherVector);
		io.write(str + "/other", otherVector, mode);
		io.write(str + "/jmPair", jmPair, mode);
		io.write(str + "/flavors", flavors, mode);
	}

	bool operator==(const Qn& a) const
	{
		return (compare(a.other) && a.oddElectrons == oddElectrons
#ifndef ENABLE_SU2
		);
#else
		    && pairEqual(a.jmPair) && flavors == a.flavors);
#endif
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

		if (direction == ProgramGlobals::DirectionEnum::INFINITE) {
			double ts = totalSites;
			for (SizeType x = 0; x < mode; ++x) {
				double flp = original.other[x] * sites;
				other[x] = static_cast<SizeType>(round(flp / ts));
			}

			double flp = original.jmPair.first * sites;
			jmPair.first = static_cast<SizeType>(round(flp / ts));
		}

		if (ifPresentOther0IsElectrons && other.size() > 0)
			oddElectrons = (other[0] & 1);

		if (!isSu2)
			return;

		assert(ifPresentOther0IsElectrons && other.size() > 0);

		SizeType tmp = jmPair.first;
		PsimagLite::String str("SymmetryElectronsSz: FATAL: Impossible parameters ");
		bool flag = false;
		if (original.oddElectrons) {
			if (!(tmp & 1)) {
				flag = true;
				str += "oddElectrons= " + ttos(original.oddElectrons) + " is odd ";
				str += "and 2j= " + ttos(tmp) + " is even.";
				tmp++;
			}
		} else {
			if (tmp & 1) {
				flag = true;
				str += "oddElectrons= " + ttos(original.oddElectrons) + " is even ";
				str += "and 2j= " + ttos(tmp) + " is odd.";
				tmp++;
			}
		}

		if (flag && sites == totalSites)
			throw PsimagLite::RuntimeError(str);

		jmPair.first = tmp;
	}

	template <typename SomeIoInType>
	static void readVector(VectorQnType& vqns,
	    PsimagLite::String prefix,
	    SomeIoInType& io)
	{
		SizeType aSize = 0;
		io.read(aSize, prefix + "/Size");

		vqns.resize(aSize, zero());
		for (SizeType i = 0; i < aSize; ++i)
			vqns[i].read(prefix + "/" + ttos(i), io);
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
			assert(1 + i * modePlusOne < ints.size());
			SizeType tmp = ints[1 + i * modePlusOne];
			assert(outQns[i].other.size() > 0);
			outQns[i].other[0] = tmp;
			outQns[i].oddElectrons = (tmp & 1);
			for (SizeType j = 1; j < modePlusOne; ++j) {
				SizeType k = (j == 1) ? 0 : j;
				assert(k + i * modePlusOne < ints.size());
				assert(j < outQns[i].other.size());
				outQns[i].other[j] = ints[k + i * modePlusOne];
			}
		}
	}

	bool isDefinedOther() const
	{
		SizeType n = other.size();
		SizeType value = 1;
		SizeType total = sizeof(value) * 8 - 1;
		value <<= total;
		for (SizeType i = 0; i < n; ++i)
			if (other[i] == value)
				return false;

		return true;
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
		SizeType total = sizeof(value) * 8 - 1;
		value <<= total;
		return Qn(false, VectorSizeType(modalStruct.size(), value), PairSizeType(0, 0), 0);
	}

	friend std::ostream& operator<<(std::ostream& os, const Qn& qn)
	{
		os << "oddElectrons=" << qn.oddElectrons << " ";
		os << "other=";
		for (SizeType i = 0; i < qn.other.size(); ++i)
			os << qn.other[i] << ",";
		os << "  jmPair=" << qn.jmPair;
		return os;
	}

	static VectorModalStructType modalStruct;
	static bool ifPresentOther0IsElectrons;
	bool oddElectrons;
	Array<SizeType> other;
	PairSizeType jmPair;
	SizeType flavors;

private:

	// assumes modulo already applied as needed
	bool compare(const Array<SizeType>& otherOther) const
	{
		SizeType n = otherOther.size();
		runChecks(n);

		for (SizeType i = 0; i < n; ++i)
			if (otherOther[i] != other[i])
				return false;

		return true;
	}

	void runChecks(SizeType n) const
	{
		assert(n == other.size());
		assert(n == modalStruct.size());
		assert(isDefinedOther());
	}

	bool pairEqual(const PairSizeType& otherJm) const
	{
		return (otherJm.first == jmPair.first && otherJm.second == jmPair.second);
	}

	void modularize()
	{
		SizeType n = other.size();
		for (SizeType i = 0; i < n; ++i) {
			if (modalStruct[i].modalEnum == MODAL_MODULO)
				other[i] %= modalStruct[i].extra;
		}
	}

	// disable implicit conversion for 1st argument of ctor
	Qn(SizeType, VectorSizeType, PairSizeType, SizeType);
};
}
#endif // QN_H
