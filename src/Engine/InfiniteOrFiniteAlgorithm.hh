#ifndef INFINITEORFINITEALGORITHM_HH
#define INFINITEORFINITEALGORITHM_HH
#include "FiniteLoop.h"
#include "ProgramGlobals.h"
#include <vector>

namespace Dmrg {
template <typename RealType> class InfiniteOrFiniteAlgorithm {
public:

	using FiniteLoopType = FiniteLoop<RealType>;

	InfiniteOrFiniteAlgorithm(ProgramGlobals::DirectionEnum      dir,
	                          const std::vector<FiniteLoopType>& finite_loops,
	                          SizeType                           loop_index)
	    : dir_(dir)
	{
		if (dir != ProgramGlobals::DirectionEnum::INFINITE) {
			if (loop_index >= finite_loops.size()) {
				err("Internal error: InfiniteOrFiniteAlgorithm\n");
			}

			this->loop_index = loop_index;
			finite_loop_     = &(finite_loops[loop_index]);
		}
	}

	bool isFinite() const { return (dir_ != ProgramGlobals::DirectionEnum::INFINITE); }

	const FiniteLoopType& finiteLoop() const
	{
		if (!finite_loop_) {
			err("InfiniteOrFiniteAlgorithm::finiteLoop() internal error\n");
		}

		return *finite_loop_;
	}

	int loopIndex() const { return loop_index; }

private:

	const ProgramGlobals::DirectionEnum dir_;
	const FiniteLoopType*               finite_loop_ = nullptr;
	int                                 loop_index   = -1;
};
}
#endif // INFINITEORFINITEALGORITHM_HH
