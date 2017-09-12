#ifndef BASESTACK_H
#define BASESTACK_H
#include "DiskStack.h"
#include <errno.h>

namespace Dmrg {

template<typename DataType>
class BaseStack {

	typedef DiskStack<DataType> DiskStackType;

public:

	BaseStack(bool disk)
	    : m_(!disk), diskStack_(0)
	{
		if (m_) return;
		diskStack_ = new DiskStackType(tmpFname(), tmpFname(),false,true);
	}

	BaseStack(const BaseStack& other)
	    : m_(other.m_), diskStack_(0)
	{
		if (m_) {
			stack_ = other.stack_;
			return;
		}

		if (other.diskStack_ == 0)
			err("other.diskStack_ pointer is 0\n");

		diskStack_ = new DiskStackType(tmpFname(), tmpFname(),false,true);
		copyDiskToDisk(*diskStack_, *(other.diskStack_));
	}

	~BaseStack()
	{
		delete diskStack_;
		diskStack_ = 0;
	}

	void push(const DataType& d)
	{
		check();
		return (m_) ? stack_.push(d) : diskStack_->push(d);
	}

	void pop()
	{
		check();
		return (m_) ? stack_.pop() : diskStack_->pop();
	}

	DataType top()
	{
		check();
		return (m_) ? stack_.top() : diskStack_->top();
	}

	SizeType size() const
	{
		check();
		return (m_) ? stack_.size() : diskStack_->size();
	}

	bool inDisk() const { return !m_; }

private:

	void check() const
	{
		if (m_) return;
		if (diskStack_ != 0) return;
		err("DiskStack pointer is 0\n");
	}

	static PsimagLite::String tmpFname()
	{
		char templ[] = "baseStackXXXXXX";
		int x = mkstemp(templ);
		if (x < 0) {
			char *msg = strerror(errno);
			err("mkstemp failed: " + PsimagLite::String(msg) + "\n");
		}
		close(x);
		return PsimagLite::String(templ);
	}

	BaseStack& operator=(const BaseStack&);

	bool m_;
	typename PsimagLite::Stack<DataType>::Type stack_;
	DiskStackType* diskStack_;
};
}
#endif // BASESTACK_H
