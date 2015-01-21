#ifndef MERSENNETWISTER_H
#define MERSENNETWISTER_H

namespace PsimagLite {

class MersenneTwister {

public:

	MersenneTwister(unsigned,int,int);

	MersenneTwister(unsigned);

	void seed(unsigned);

	unsigned random();

	double operator() ();

	static unsigned max() { return keepLast32BitMask_; }

private:

	static const unsigned N_ = 624;
	static const unsigned keepLast32BitMask_ = 4294967295; // 2^32 - 1
	unsigned index_;
	unsigned state_[N_];
	void generate();
};

} // namespace PsimagLite

/*  K. A. Al-Hassanieh
 *
 *  MersenneTwister generates random 32-bit unsinged integers.
 *  The constructor takes an unsigned int as a seed.
 *  The return type of random() is unsigned int. Converting the
 *  output to int can cause errors, since the returned value can
 *  be larger than the max of signed int  which is 2^16 - 1:
 *  int i = rng.random(); // causes errors
 *  instead use:
 *  unsigned i = rng.random();
 *  or
 *  size_t i = rng.random();
 *
 *  operator() returns a double between 0 and 1.
 */

#endif // MERSENNETWISTER_H

