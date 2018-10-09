#ifndef DISKORMEMORYSTACK_H
#define DISKORMEMORYSTACK_H
#include "Stack.h"
#include "DiskStackNg.h"

namespace Dmrg {

template<typename BasisWithOperatorsType>
class DiskOrMemoryStack {

public:

	typedef typename PsimagLite::Stack<BasisWithOperatorsType>::Type MemoryStackType;

	void push(const BasisWithOperatorsType& b) { memory_.push(b); }

	void pop() { memory_.pop(); }

	SizeType size() const { return memory_.size(); }

	const BasisWithOperatorsType& top() const { return memory_.top(); }

private:

	MemoryStackType memory_;
};
}
#endif // DISKORMEMORYSTACK_H
