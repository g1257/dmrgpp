#ifndef DISKORMEMORYSTACK_H
#define DISKORMEMORYSTACK_H
#include "Stack.h"
#include "DiskStackNg.h"
#include "Io/IoNg.h"

namespace Dmrg {

template<typename BasisWithOperatorsType>
class DiskOrMemoryStack {

public:

	typedef typename PsimagLite::Stack<BasisWithOperatorsType>::Type MemoryStackType;
	typedef DiskStack<BasisWithOperatorsType> DiskStackType;

	DiskOrMemoryStack(bool onDisk,
	                  const PsimagLite::String filename,
	                  PsimagLite::String label,
	                  bool isObserveCode)
	    : diskW_(0), diskR_(0)
	{
		if (!onDisk) return;

		size_t lastindex = filename.find_last_of(".");
		PsimagLite::String file = filename.substr(0, lastindex) + "Stacks.hd5";

		if (createFile_) {
			PsimagLite::IoNg::Out out(file, PsimagLite::IoNg::ACC_TRUNC);
			out.close();
			createFile_ = false;
		}

		diskW_ = new DiskStackType(file, false, label, isObserveCode);
		diskR_ = new DiskStackType(file, true, label, isObserveCode);
	}

	~DiskOrMemoryStack()
	{
		delete diskR_;
		diskR_ = 0;
		delete diskW_;
		diskW_ = 0;
	}

	void push(const BasisWithOperatorsType& b)
	{
		if (diskW_) {
			diskW_->push(b);
			diskW_->flush();
			diskR_->restore(diskW_->size());
		} else {
			memory_.push(b);
		}
	}

	void pop()
	{
		if (diskW_) {
			diskW_->pop();
			diskW_->flush();
			diskR_->restore(diskW_->size());
		} else {
			memory_.pop();
		}
	}

	bool onDisk() const { return (diskR_); }

	SizeType size() const
	{
		return (diskR_) ? diskR_->size() : memory_.size();
	}

	const BasisWithOperatorsType& top() const
	{
		return (diskR_) ? diskR_->top() : memory_.top();
	}

	void toDisk(DiskStackType& disk) const
	{
		if (diskR_) {
			SizeType total = diskR_->size();
			DiskStackType& diskNonConst = const_cast<DiskStackType&>(*diskR_);
			loadStack(disk, diskNonConst);
			assert(diskW_);
			diskW_->restore(total);
		} else {
			MemoryStackType memory = memory_;
			loadStack(disk, memory);
		}
	}

	template<typename StackType1,typename StackType2>
	static void loadStack(StackType1& stackInMemory, StackType2& stackInDisk)
	{
		while (stackInDisk.size()>0) {
			BasisWithOperatorsType b = stackInDisk.top();
			stackInMemory.push(b);
			stackInDisk.pop();
		}
	}

	static bool createFile_;

private:

	DiskOrMemoryStack(const DiskOrMemoryStack&);

	DiskOrMemoryStack& operator=(const DiskOrMemoryStack&);

	MemoryStackType memory_;
	DiskStackType *diskW_;
	DiskStackType *diskR_;
};

template<typename BasisWithOperatorsType>
bool DiskOrMemoryStack<BasisWithOperatorsType>::createFile_ = true;
}
#endif // DISKORMEMORYSTACK_H
