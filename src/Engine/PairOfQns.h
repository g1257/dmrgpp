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

	Qn make() const
	{
		assert(q1_ && q2_);
		return Qn(*q1_, *q2_);
	}

private:

	const QnType* q1_;
	const QnType* q2_;
};
}
#endif // PAIROFQNS_H
