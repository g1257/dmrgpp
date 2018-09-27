#ifndef PAIROFQNS_H
#define PAIROFQNS_H
#include "Qn.h"
namespace Dmrg {

class PairOfQns {

public:

	typedef Qn QnType;

	PairOfQns() : q1_(0), q2_(0) {}

	PairOfQns(const QnType& q1, const QnType& q2)
	    : q1_(&q1), q2_(&q2)
	{}

	bool operator==(const PairOfQns& other) const
	{
		assert(q1_ && q2_);
		const QnType& q1 = *q1_;
		const QnType& q2 = *q2_;
		SizeType n = q1.other.size();
		assert(q2.other.size() == n);
		assert(Qn::modalStruct.size() == n);

		assert(other.q1_ && other.q2_);
		const QnType& r1 = *other.q1_;
		const QnType& r2 = *other.q2_;
		assert(r1.other.size() == n);
		assert(r2.other.size() == n);

		for (SizeType i = 0; i < n; ++i) {
			SizeType a = q1.other[i] + q2.other[i];
			SizeType b = r1.other[i] + r2.other[i];

			if (Qn::modalStruct[i].modalEnum == Qn::MODAL_MODULO) {
				a %= Qn::modalStruct[i].extra;
				b %= Qn::modalStruct[i].extra;
			}

			if (a != b) return false;
		}

		bool oa = (q1.oddElectrons ^ q2.oddElectrons);
		bool ob = (r1.oddElectrons ^ r2.oddElectrons);

#if ENABLE_SU2
		assert(false);
#endif

		return (oa == ob);
	}

	SizeType hash(bool addOdd) const
	{
		assert(q1_ && q2_);
		const QnType& q1 = *q1_;
		const QnType& q2 = *q2_;
		bool oddElectrons = (q1.oddElectrons ^ q2.oddElectrons);
		SizeType n = q1.other.size();
		assert(q2.other.size() == n);
		assert(Qn::modalStruct.size() == n);

		const SizeType offset = 8; // 8 bits
		SizeType key = (addOdd && oddElectrons) ? 1 : 0;
		SizeType bits = (addOdd) ? 1 : 0;

		for (SizeType i = 0; i < n; ++i) {
			SizeType val = q1.other[i] + q2.other[i];
			if (Qn::modalStruct[i].modalEnum == Qn::MODAL_MODULO)
				val %= Qn::modalStruct[i].extra;
			val <<= bits;
			key += val;
			bits += offset;
		}

		return key;
	}

	bool oddElectrons() const
	{
		assert(q1_ && q2_);
		return (q1_->oddElectrons ^ q2_->oddElectrons);
	}

	// expensive
	Qn make() const
	{
		assert(q1_ && q2_);
		return Qn(*q1_, *q2_);
	}

private:

	//PairOfQns(const PairOfQns&);

	//PairOfQns& operator=(const PairOfQns&);

	const QnType* q1_;
	const QnType* q2_;
};
}
#endif // PAIROFQNS_H
