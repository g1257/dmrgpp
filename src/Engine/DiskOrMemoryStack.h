#ifndef DISKORMEMORYSTACK_H
#define DISKORMEMORYSTACK_H
#include "DiskStackNg.h"
#include "Io/IoNg.h"
#include "Stack.h"

namespace Dmrg
{

template <typename T>
class MemoryStack : public PsimagLite::Stack<T>::Type
{

public:

	using PsimagLite::Stack<T>::Type::c;
};

template <typename BasisWithOperatorsType>
class DiskOrMemoryStack
{

public:

	typedef MemoryStack<BasisWithOperatorsType> MemoryStackType;
	typedef DiskStack<BasisWithOperatorsType> DiskStackType;

	DiskOrMemoryStack(bool onDisk,
	    const PsimagLite::String filename,
	    const PsimagLite::String post,
	    PsimagLite::String label,
	    const BasisTraits& basisTraits)
	    : basisTraits_(basisTraits)
	    , diskW_(0)
	    , diskR_(0)
	{
		if (!onDisk)
			return;

		size_t lastindex = filename.find_last_of(".");
		PsimagLite::String file = filename.substr(0, lastindex) + post + ".hd5";

		if (createFile_) {
			PsimagLite::IoNg::Out out(file, PsimagLite::IoNg::ACC_TRUNC);
			out.close();
			createFile_ = false;
		}

		diskW_ = new DiskStackType(file, false, label, basisTraits_);
		diskR_ = new DiskStackType(file, true, label, basisTraits_);
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
			diskR_->restore(total);
		} else {
			MemoryStackType memory = memory_;
			loadStack(disk, memory);
		}
	}

	void read(PsimagLite::String prefix, PsimagLite::IoNgSerializer& io)
	{
		if (diskW_) {
			assert(diskR_);
			readWftStacksOnDisk(prefix, io);
		}

		io.read(memory_, prefix);
	}

	void write(PsimagLite::String prefix, PsimagLite::IoNgSerializer& io) const
	{
		if (diskW_) {
			assert(diskR_);
			const SizeType total = diskR_->size();
			writeWftStacksOnDisk(prefix, io);
			diskR_->restore(total);
			diskW_->restore(total);
			return;
		}

		MemoryStackType m = memory_;
		io.write(prefix, m);
	}

	const BasisWithOperatorsType& operator[](SizeType ind) const
	{
		if (diskW_)
			err("MultipointInSitu does not support stacks on disk, only on memory\n");

		return memory_.c[ind];
	}

	template <typename StackType1, typename StackType2>
	static void loadStack(StackType1& stackInMemory, StackType2& stackInDisk)
	{
		while (stackInDisk.size() > 0) {
			BasisWithOperatorsType b = stackInDisk.top();
			stackInMemory.push(b);
			stackInDisk.pop();
		}
	}

	static bool createFile_;

private:

	void writeWftStacksOnDisk(PsimagLite::String name,
	    PsimagLite::IoNgSerializer& io) const
	{
		io.createGroup(name);
		io.write(name + "/Size", this->size());
		SizeType i = 0;
		while (this->size() > 0) {
			const BasisWithOperatorsType& t = this->top();
			t.write(name + "/" + ttos(i++), io);
			DiskOrMemoryStack& thisObject = const_cast<DiskOrMemoryStack&>(*this);
			thisObject.pop();
		}
	}

	void readWftStacksOnDisk(PsimagLite::String name, PsimagLite::IoNgSerializer& io)
	{
		SizeType x = 0;
		io.read(x, name + "/Size");
		for (SizeType i = 0; i < x; ++i) {
			BasisWithOperatorsType t;
			t.read(name + "/" + ttos(x - i - 1), io);
			this->push(t);
		}
	}

	DiskOrMemoryStack(const DiskOrMemoryStack&);

	DiskOrMemoryStack& operator=(const DiskOrMemoryStack&);

	const BasisTraits& basisTraits_;
	mutable MemoryStackType memory_;
	DiskStackType* diskW_;
	DiskStackType* diskR_;
};

template <typename BasisWithOperatorsType>
bool DiskOrMemoryStack<BasisWithOperatorsType>::createFile_ = true;
}
#endif // DISKORMEMORYSTACK_H
