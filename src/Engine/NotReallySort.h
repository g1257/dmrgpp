#ifndef NOT_REALLY_SORT_H
#define NOT_REALLY_SORT_H
#define USE_PTHREADS_OR_NOT_NG
#include "Qn.h"
#include "Vector.h"
#include <numeric>
#include "Sort.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include <unordered_map>
#include "PairOfQns.h"
#include "Array.h"

namespace std {

template<>
class hash<Dmrg::Qn> {

public:

	typedef Dmrg::Qn::VectorQnType VectorQnType;
	typedef Dmrg::Qn::VectorSizeType VectorSizeType;

	hash(VectorSizeType* hash,
	     const VectorQnType* inQns,
	     bool addOdd)
	    : hash_(hash), inQns_(inQns), addOdd_(addOdd)
	{}

	hash(bool addOdd) : addOdd_(addOdd)
	{}

	void doTask(SizeType taskNumber, SizeType)
	{
		assert(inQns_->size() > taskNumber);
		(*hash_)[taskNumber] = operator()((*inQns_)[taskNumber]);
	}

	SizeType tasks() const
	{
		return inQns_->size();
	}

	SizeType operator()(const Dmrg::Qn& qn) const
	{
		const SizeType offset = 8; // 8 bits
		const SizeType n = qn.other.size(); // small number
		SizeType key = (addOdd_ && qn.oddElectrons) ? 1 : 0;
		SizeType bits = (addOdd_) ? 1 : 0;

		for (SizeType i = 0; i < n; ++i) {
			SizeType val = qn.other[i];
			val <<= bits;
			key += val;
			bits += offset;
		}

		return key;
	}

private:

	VectorSizeType* hash_;
	const VectorQnType* inQns_;
	bool addOdd_;
};

template<>
class hash<Dmrg::PairOfQns> {

public:

	typedef Dmrg::Array<Dmrg::PairOfQns> VectorLikeQnType;
	typedef Dmrg::Qn::VectorSizeType VectorSizeType;


	hash(VectorSizeType& hash,
	     const VectorLikeQnType& inQns,
	     bool addOdd)
	    : hash_(hash), inQns_(inQns), addOdd_(addOdd)
	{}

	SizeType operator()(const Dmrg::PairOfQns& qnPair) const
	{
		return qnPair.hash(addOdd_);
	}

	void doTask(SizeType taskNumber, SizeType)
	{
		assert(inQns_.size() > taskNumber);
		hash_[taskNumber] = operator()(inQns_[taskNumber]);
	}

	SizeType tasks() const
	{
		return inQns_.size();
	}

private:

	VectorSizeType& hash_;
	const VectorLikeQnType& inQns_;
	bool addOdd_;
};
} // namespace std

#endif // NOT_REALLY_SORT_H

