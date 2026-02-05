#ifndef DMRG_QN_HASH_H
#define DMRG_QN_HASH_H
#define USE_PTHREADS_OR_NOT_NG
#include "Array.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "Qn.h"
#include "Sort.h"
#include "Vector.h"
#include <numeric>
#include <unordered_map>

namespace std {

template <> class hash<Dmrg::Qn> {

public:

	using VectorQnType   = Dmrg::Qn::VectorQnType;
	using VectorSizeType = Dmrg::Qn::VectorSizeType;

	hash(bool addOdd)
	    : addOdd_(addOdd)
	{ }

	SizeType operator()(const Dmrg::Qn& qn) const
	{
		const SizeType offset = 8; // 8 bits
		const SizeType n      = qn.other.size(); // small number
		SizeType       key    = (addOdd_ && qn.oddElectrons) ? 1 : 0;
		SizeType       bits   = (addOdd_) ? 1 : 0;

		for (SizeType i = 0; i < n; ++i) {
			SizeType val = qn.other[i];
			val <<= bits;
			key += val;
			bits += offset;
		}

		return key;
	}

private:

	bool addOdd_;
};

} // namespace std

#endif // DMRG_QN_HASH_H
