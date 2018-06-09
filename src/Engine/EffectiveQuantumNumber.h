#ifndef EFFECTIVEQUANTUMNUMBER_H
#define EFFECTIVEQUANTUMNUMBER_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Io/IoSelector.h"
#include "TargetQuantumElectrons.h"
#include "Io/IoNgSerializer.h"
#include "Sort.h"
#include "SymmetryElectronsSz.h"

namespace Dmrg {

template<typename RealType>
class EffectiveQuantumNumber {

	typedef PsimagLite::IoSelector::Out IoOutType;

	class Qn {

	public:

		void read(PsimagLite::String str, PsimagLite::IoNg::In& io)
		{
			io.read(data_, str);
		}

		void read(PsimagLite::String str, PsimagLite::IoNgSerializer& io)
		{
			io.read(data_, str);
		}

		void write(PsimagLite::String str, PsimagLite::IoNgSerializer& io) const
		{
			io.write(str, data_);
		}

		bool operator==(const Qn& a) const
		{
			return (a.data_ == data_);
		}

		bool operator!=(const Qn& a) const
		{
			return !(*this == a);
		}

		SizeType toInteger() const
		{
			return data_;
		}

		friend class EffectiveQuantumNumber;

		friend std::ostream& operator<<(std::ostream& os, const Qn& qn)
		{
			os<<qn.data_;
			return os;
		}

	private:

		SizeType data_;
	};

public:

	typedef Qn QnType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<QnType>::Type VectorQnType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;
	typedef std::pair<SizeType, SizeType> PairType;
	typedef SymmetryElectronsSz<RealType> SymmetryElectronsSzType;

	static QnType tensorProduct(const QnType& q1, const QnType& q2)
	{
		QnType q;
		q.data_ = q1.data_ + q2.data_;
		return q;
	}

	static QnType fromInteger(SizeType x)
	{
		QnType q;
		q.data_ = x;
		return q;
	}

	static void findQuantumNumbers(VectorSizeType& qn,
	                               const SymmetryElectronsSzType& qq,
	                               bool useSu2Symmetry)
	{
		if (useSu2Symmetry)
			findQuantumNumbersSu2(qn, qq);
		else
			findQuantumNumbersLocal(qn, qq);
	}

	static void sort(VectorSizeType& q, VectorSizeType& iperm)
	{
		PsimagLite::Sort<VectorSizeType> sort;
		sort.sort(q, iperm);
	}

	static PsimagLite::String qnPrint(const QnType& q, SizeType total)
	{
		PsimagLite::String str("");
		VectorSizeType qns = decodeQuantumNumber(q.data_, total);
		for (SizeType k = 0; k < qns.size(); ++k)
			str += ttos(qns[k]) + " ";
		return str;
	}

	static QnType adjustQn(const VectorSizeType& adjustQuantumNumbers,
	                       ProgramGlobals::DirectionEnum direction,
	                       IoOutType& ioOut,
	                       bool useSu2Symmetry,
	                       SizeType step,
	                       SizeType mode)
	{
		VectorSizeType targetQuantumNumbers(mode+1,0);
		if (2*step+mode >= adjustQuantumNumbers.size()) {
			PsimagLite::String msg("adjustQuantumNumbers must be a vector");
			msg += " of correct size\n";
			throw PsimagLite::RuntimeError(msg);
		}

		for (SizeType x = 0; x < (mode+1); ++x)
			targetQuantumNumbers[x] = adjustQuantumNumbers[2*step+x];

		return getQuantumSector(targetQuantumNumbers,
		                        direction,
		                        &ioOut,
		                        useSu2Symmetry);
	}

	static QnType getQuantumSector(const TargetQuantumElectronsType& targetQuantum,
	                               SizeType sites,
	                               SizeType total,
	                               ProgramGlobals::DirectionEnum direction,
	                               IoOutType* ioOut,
	                               bool useSu2Symmetry)
	{
		VectorSizeType v;
		setTargetNumbers(v,targetQuantum,sites,total,direction);
		return getQuantumSector(v,direction,ioOut,useSu2Symmetry);
	}

	static void qnToElectrons(VectorSizeType& electrons,
	                          const VectorSizeType& qns,
	                          SizeType total)
	{
		electrons.resize(qns.size());
		for (SizeType i=0;i<qns.size();i++) {
			VectorSizeType v = decodeQuantumNumber(qns[i],total);
			electrons[i] = v[1];
		}
	}

	static QnType pseudoEffectiveNumber(SizeType nelectrons,
	                                      SizeType jtilde)
	{
		VectorSizeType v(3);
		v[0] = 0;
		v[1] = nelectrons;
		v[2] = jtilde;
		return fromInteger(pseudoQuantumNumber_(v));
	}

	static SizeType neJmToIndex(SizeType ne,const PairType& jm)
	{
		VectorSizeType v(3);
		v[0] = jm.second;
		v[1] = ne;
		v[2] = jm.first;
		return encodeQuantumNumber(v);
	}

private:

	static void findQuantumNumbersSu2(VectorSizeType& q, const SymmetryElectronsSzType& qq)
	{
		q.resize(qq.electrons.size());
		for (SizeType i=0;i<q.size();i++) {
			SizeType ne = qq.electrons[i];
			PairType jmpair = qq.jmValues[i];
			q[i]=neJmToIndex(ne,jmpair);
		}
	}

	//! find quantum numbers for each state of this basis,
	//! considered symmetries for this model are: n_up and n_down
	static void findQuantumNumbersLocal(VectorSizeType& q, const SymmetryElectronsSzType& qq)
	{
		SizeType mode = static_cast<SizeType>(qq.other.size()/qq.electrons.size());
		assert(qq.other.size() % qq.electrons.size() == 0);
		assert(mode > 0);

		q.clear();
		q.resize(qq.electrons.size());

		VectorSizeType qn(mode + 1);
		for (SizeType i = 0 ; i < qq.electrons.size(); ++i) {
			// n
			qn[1] = qq.electrons[i];
			// sz + const.
			qn[0] = qq.other[i];

			for (SizeType x = 0; x < (mode-1); ++x)
				qn[2+x] = qq.other[i + (x+1)*qq.electrons.size()];
			q[i] = encodeQuantumNumber(qn);
		}
	}

	static void setTargetNumbers(VectorSizeType& t,
	                             const TargetQuantumElectronsType& targetQ,
	                             SizeType sites,
	                             SizeType totalSites,
	                             ProgramGlobals::DirectionEnum direction)
	{
		SizeType mode = targetQ.other.size();
		assert(mode > 0);
		assert(!targetQ.isSu2 || mode == 1);

		t.resize((targetQ.isSu2) ? 3 : mode + 1,0);

		if (direction == ProgramGlobals::INFINITE) {
			t[0] = static_cast<SizeType>(round(targetQ.other[0]*sites/totalSites));
			t[1] = static_cast<SizeType>(round(targetQ.totalElectrons*sites/totalSites));
			for (SizeType x = 0; x < (mode-1); ++x)
				t[2+x] = static_cast<SizeType>(round(targetQ.other[x+1]*sites/totalSites));
		} else {
			t[0] = targetQ.other[0];
			t[1] = targetQ.totalElectrons;
			for (SizeType x = 0; x < (mode-1); ++x)
				t[2+x] = static_cast<SizeType>(round(targetQ.other[x+1]*sites/totalSites));
		}

		if (!targetQ.isSu2) return;

		RealType jReal = targetQ.twiceJ*sites/static_cast<RealType>(totalSites);
		SizeType tmp = (direction == ProgramGlobals::INFINITE) ?
		            static_cast<SizeType>(round(jReal)) : targetQ.twiceJ;

		PsimagLite::String str("SymmetryElectronsSz: FATAL: Impossible parameters ");
		bool flag = false;
		if (targetQ.totalElectrons & 1) {
			if (!(tmp&1)) {
				flag = true;
				str += "electrons= " + ttos(targetQ.totalElectrons) + " is odd ";
				str += "and 2j= " +  ttos(tmp) + " is even.";
				tmp++;
			}
		} else {
			if (tmp & 1) {
				flag = true;
				str += "electrons= " + ttos(targetQ.totalElectrons) + " is even ";
				str += "and 2j= " +  ttos(tmp) + " is odd.";
				tmp++;
			}
		}

		if (flag && sites == totalSites) throw PsimagLite::RuntimeError(str);

		t[2] = tmp;
	}

	static VectorSizeType decodeQuantumNumber(SizeType q, SizeType total)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		SizeType number = 1;
		for (SizeType x = 1; x < total; ++x) number *= maxElectrons;
		SizeType tmp = q;
		VectorSizeType v(total);
		for (SizeType x = 0; x < total; ++x) {
			v[total-1-x] = static_cast<SizeType>(tmp/number);
			tmp -= v[total-1-x]*number;
			number /= maxElectrons;

		}

		return v;
	}

	static QnType getQuantumSector(const VectorSizeType& targetQuantumNumbers,
	                               ProgramGlobals::DirectionEnum direction,
	                               IoOutType* ioOut,
	                               bool useSu2Symmetry)
	{
		static SizeType counter = 0;
		PsimagLite::OstringStream msg;
		msg<<"SymmetryElectronsSz: Integer target quantum numbers are: ";
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		std::cout<<msg.str()<<"\n";
		if (ioOut && direction == ProgramGlobals::INFINITE)
			ioOut->write(targetQuantumNumbers,"TargetedQuantumNumbers" + ttos(counter));
		++counter;
		return fromInteger(pseudoQuantumNumber(targetQuantumNumbers,useSu2Symmetry));
	}

	//! Encodes (flavor,jvalue,density) into a unique number and returns it
	static SizeType pseudoQuantumNumber(const VectorSizeType& targets,
	                                    bool useSu2Symmetry)
	{
		if (useSu2Symmetry)
			return pseudoQuantumNumber_(targets);
		else
			return encodeQuantumNumber(targets);
	}

	//! targets[0]=nup, targets[1]=ndown,  targets[2]=2j
	static SizeType pseudoQuantumNumber_(const VectorSizeType& v)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		SizeType x = v[1];

		assert(x < maxElectrons);

		x += v[2]*maxElectrons;
		return x;
	}

	static SizeType encodeQuantumNumber(const VectorSizeType& v)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		assert(v.size() > 0);
		SizeType n = (v.size() == 0) ? 0 : v.size() - 1;
		for (SizeType x = 0; x < n; ++x)
			assert(v[x] < maxElectrons);

		SizeType index = 0;
		SizeType number = 1;
		for (SizeType x = 0; x < v.size(); ++x) {
			index += v[x] * number;
			number *= maxElectrons;
		}

		return index;
	}

	PsimagLite::String str_;
};
}
#endif // EFFECTIVEQUANTUMNUMBER_H
