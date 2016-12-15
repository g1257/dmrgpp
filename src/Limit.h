#ifndef PSIMAGLITE_LIMIT_H
#define PSIMAGLITE_LIMIT_H
#include <sys/time.h>
#include <sys/resource.h>
#include "Vector.h"
#include <utility>
#include <string.h>
#include <errno.h>

namespace PsimagLite {

class Limit {

public:

	typedef std::pair<rlim_t, rlim_t> PairRlimType;

	Limit() : rlimit_(new struct rlimit)
	{}

	~Limit()
	{
		delete rlimit_;
		rlimit_ = 0;
	}

	void memory(int s, int h = 0)
	{
		rlimit_->rlim_cur = s;
		rlimit_->rlim_max = (h > s) ? h : s;
		int ret = setrlimit(RLIMIT_AS,rlimit_);
		checkRet(ret,"setrlimit");
	}

	PairRlimType memory()
	{
		int ret = getrlimit(RLIMIT_AS,rlimit_);
		checkRet(ret,"getrlimit");
		return PairRlimType(rlimit_->rlim_cur,rlimit_->rlim_max);
	}

private:

	void checkRet(int x, PsimagLite::String msg) const
	{
		if (x == 0) return;
		std::cerr<<"Call to "<<msg<<" failed\n";
		std::cerr<<strerror(errno)<<"\n";
		// FIXME: should we throw here?
	}

	struct rlimit* rlimit_;
};
}; // namespace PsimagLite
#endif

