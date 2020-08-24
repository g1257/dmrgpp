#ifndef OUTPUTFILEORNOT_H
#define OUTPUTFILEORNOT_H
#include "Io/IoNg.h"
#include "Io/IoSelector.h"

namespace Dmrg {

class OutputFileOrNot {

public:

	OutputFileOrNot(PsimagLite::String filename,
	                PsimagLite::IoNg::OpenMode mode,
	                bool disabled)
	    : filename_(filename), ptr_(nullptr)
	{
		if (!disabled)
			ptr_ = new PsimagLite::IoSelector::Out(filename, mode);
	}

	~OutputFileOrNot()
	{
		delete ptr_;
		ptr_ = 0;
	}

	const PsimagLite::String& filename() const { return filename_; }

	bool nonNull() const { return (ptr_ != nullptr); }

	PsimagLite::IoSelector::Out& handle()
	{
		if (!ptr_)
			err("OutputFileOrNot: FATAL: Requesting handle of empty object\n");
		return *ptr_;
	}

	template<typename T>
	void write(const T& t, PsimagLite::String str)
	{
		if (!ptr_) return;
		ptr_->write(t, str);
	}

	void write(SizeType c,
	           PsimagLite::String str,
	           PsimagLite::IoNg::Out::Serializer::WriteMode mode)
	{
		if (!ptr_) return;
		ptr_->write(c, str, mode);
	}

	void flush()
	{
		if (!ptr_) return;
		ptr_->flush();
	}

	void close()
	{
		if (!ptr_) return;
		ptr_->close();
	}

	void createGroup(PsimagLite::String str)
	{
		if (!ptr_) return;
		ptr_->createGroup(str);
	}

private:

	PsimagLite::String filename_;
	PsimagLite::IoSelector::Out* ptr_;
};
}
#endif // OUTPUTFILEORNOT_H
