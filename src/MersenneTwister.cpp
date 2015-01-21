#include "MersenneTwister.h"

namespace PsimagLite {

MersenneTwister::MersenneTwister(unsigned s)
    : index_(0)
{
	seed(s);
}

MersenneTwister::MersenneTwister(unsigned s, int rank, int nprocs)
: index_(0)
{
	seed(s);
	int news = random();
	for (int i = 1; i < rank; ++i) {
		news = random();
	}

	seed(news);
}

void
MersenneTwister::seed(unsigned s)
{
	index_ = 0;
	state_[0] = s;
	for (unsigned i = 1; i < N_; ++i) {
		unsigned tmp = state_[i-1]^(state_[i-1]>>30);	
		state_[i] =  keepLast32BitMask_ &
		        (1812433253 * (tmp + i));
	}
}

void
MersenneTwister::generate()
{
	for (unsigned i = 0; i < N_; ++i) {
		unsigned y = (state_[i] & 0x80000000) +
		        (state_[(i+1)%N_] & 0x7fffffff);
		state_[i] = state_[(i+397)%N_] ^ (y >> 1);
		if (y%2 != 0)
			state_[i] ^= 2567483615;
	}
}

unsigned
MersenneTwister::random()
{
	if (index_ == 0)
		generate();

	unsigned y = state_[index_];
	y ^= (y >> 11);
	y ^= (y << 07) & 2636928640;
	y ^= (y << 15) & 4022730752;
	y ^= (y >> 18);

	index_ = (index_+1)%N_;

	return y;
}

double
MersenneTwister::operator() ()
{
	unsigned r = random();
	return (static_cast<double>(r)/keepLast32BitMask_);
}

} // namespace PsimagLite

