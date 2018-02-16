#ifndef BASESTACK_H
#define BASESTACK_H
#include "DiskStack.h"
#include <errno.h>
#include "IoSelector.h"

namespace Dmrg {

template<typename DataType>
class BaseStack {

	typedef DiskStack<DataType> DiskStackType;
	typedef PsimagLite::Vector<PsimagLite::String>::Type VectorStringType;
	typedef typename PsimagLite::IoSelector::In IoInType;
	typedef typename PsimagLite::IoSelector::Out IoOutType;

public:

	BaseStack(bool disk)
	    : m_(!disk), diskStack_(0)
	{
		if (m_) return;
		PsimagLite::String tmpfname = tmpFname();
		files_.push_back(tmpfname);
		diskStack_ = new DiskStackType(tmpfname, tmpfname, false, false);
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

		PsimagLite::String tmpfname = tmpFname();
		files_.push_back(tmpfname);
		diskStack_ = new DiskStackType(tmpfname, tmpfname, false, false);
		copyDiskToDisk(*diskStack_, *(other.diskStack_));
	}

	~BaseStack()
	{
		deleteFiles();

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

	const DataType& top() const
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

	void save(IoOutType& io, PsimagLite::String label) const
	{
		if (m_) {
			io.print(label, stack_);
			return;
		}

		diskStack_->copyToIo(io, label);
	}

	void load(IoInType& io, PsimagLite::String label)
	{
		if (m_) {
			io.read(stack_, label);
			return;
		}

		diskStack_->copyFromIo(io, label);
	}

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

	void deleteFiles()
	{
		for (SizeType i = 0; i < files_.size(); ++i) {
			int x = unlink(files_[i].c_str());
			if (x == 0) continue;
			std::cerr<<"unlink "<<files_[i]<<" failed\n";
			std::cerr<<strerror(errno)<<"\n";
		}
	}

	BaseStack& operator=(const BaseStack&);

	bool m_;
	typename PsimagLite::Stack<DataType>::Type stack_;
	DiskStackType* diskStack_;
	VectorStringType files_;
};
}
#endif // BASESTACK_H
