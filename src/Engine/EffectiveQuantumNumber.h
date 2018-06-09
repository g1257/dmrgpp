#ifndef EFFECTIVEQUANTUMNUMBER_H
#define EFFECTIVEQUANTUMNUMBER_H
#include "Vector.h"
#include "ProgramGlobals.h"
#include "Io/IoSelector.h"
#include "TargetQuantumElectrons.h"
#include "Io/IoNgSerializer.h"

namespace Dmrg {

template<typename RealType_>
class EffectiveQuantumNumber {

	typedef PsimagLite::IoSelector::Out IoOutType;

	class Qn {

	public:

		void read(PsimagLite::String str, PsimagLite::IoNg::In& io)
		{
			throw PsimagLite::RuntimeError("unimplemented\n");
		}

		void read(PsimagLite::String str, PsimagLite::IoNgSerializer& io)
		{
			throw PsimagLite::RuntimeError("unimplemented\n");
		}

		void write(PsimagLite::String str, PsimagLite::IoNgSerializer& io) const
		{
			throw PsimagLite::RuntimeError("unimplemented\n");
		}

		bool operator==(const Qn& a) const
		{
			err("Unimplemented\n");
			return false;
		}

		bool operator!=(const Qn& a) const
		{
			return !(*this == a);
		}

		bool isDefined() const
		{
			throw PsimagLite::RuntimeError("unimplemented\n");
		}

		SizeType toInteger() const
		{
			throw PsimagLite::RuntimeError("unimplemented\n");
		}

		friend std::ostream& operator<<(std::ostream& os, const Qn& qn)
		{
			throw PsimagLite::RuntimeError("unimplemented\n");
		}

	private:

		unsigned char data_[8];
	};
public:

	typedef RealType_ RealType;
	typedef Qn QnType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<QnType>::Type VectorQnType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;
	typedef std::pair<SizeType, SizeType> PairType;

	const PsimagLite::String& name() const { return str_;}

	static QnType tensorProduct(const QnType& q1, const QnType& q2)
	{
		throw PsimagLite::RuntimeError("unimplemented\n");
	}

	static QnType fromInteger(SizeType)
	{
		throw PsimagLite::RuntimeError("unimplemented\n");
	}

	static void sort(VectorSizeType& q, VectorSizeType& iperm)
	{
		err("Unimplemented\n");
	}

	template<typename T>
	static void findQuantumNumbers(VectorSizeType& q, T& s, bool)
	{
		err("Unimplemented\n");
	}

	static PsimagLite::String qnPrint(const QnType& q, SizeType total)
	{
		return "unimplemented";
		//		PsimagLite::String str("");
		//		VectorSizeType qns = decodeQuantumNumber(q,total);
		//		for (SizeType k=0;k<qns.size();k++) str += ttos(qns[k]) + " ";
		//		return str;
	}

	static QnType adjustQn(const VectorSizeType& adjustQuantumNumbers,
	                       ProgramGlobals::DirectionEnum direction,
	                       IoOutType& ioOut,
	                       bool useSu2Symmetry,
	                       SizeType step,
	                       SizeType mode)
	{
		throw PsimagLite::RuntimeError("unimplemented\n");
	}

	static QnType getQuantumSector(const TargetQuantumElectronsType& targetQuantum,
	                               SizeType sites,
	                               SizeType total,
	                               ProgramGlobals::DirectionEnum direction,
	                               IoOutType* ioOut,
	                               bool useSu2Symmetry)
	{
		throw PsimagLite::RuntimeError("unimplemented\n");
	}

	static void qnToElectrons(VectorSizeType& electrons,
	                          const VectorSizeType& qns,
	                          SizeType total)
	{
		//		electrons.resize(qns.size());
		//		for (SizeType i=0;i<qns.size();i++) {
		//			VectorSizeType v = decodeQuantumNumber(qns[i],total);
		//			electrons[i] = v[1];
		//		}
		throw PsimagLite::RuntimeError("unimplemented\n");
	}

	static SizeType neJmToIndex(SizeType ne,const PairType& jm)
	{
		throw PsimagLite::RuntimeError("unimplemented\n");
//		VectorSizeType v(3);
//		v[0] = jm.second;
//		v[1] = ne;
//		v[2] = jm.first;
//		return encodeQuantumNumber(v);
	}

	static QnType pseudoEffectiveNumber(SizeType nelectrons,
	                                    SizeType jtilde)
	{
		//		VectorSizeType v(3);
		//		v[0] = 0;
		//		v[1] = nelectrons;
		//		v[2] = jtilde;
		//		return pseudoQuantumNumber_(v);
		throw PsimagLite::RuntimeError("unimplemented\n");

	}

private:

	PsimagLite::String str_;
};
}
#endif // EFFECTIVEQUANTUMNUMBER_H
