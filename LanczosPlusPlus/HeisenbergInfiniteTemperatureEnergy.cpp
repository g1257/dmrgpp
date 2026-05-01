#include <cassert>
#include <iostream>

bool isValid(unsigned int x, unsigned int twiceS) { return (x < twiceS + 1); }

int findSzPlusConst(int x, int n, int oneSiteSizeInBits, unsigned int twiceS)
{
	unsigned int mask = (1 << oneSiteSizeInBits) - 1;
	int          sum  = 0;
	for (int i = 0; i < n; ++i) {
		unsigned int val = (x & mask);
		if (!isValid(val, twiceS))
			return -1;
		sum += val;
		x >>= oneSiteSizeInBits;
	}

	assert(x == 0);
	return sum;
}

double energy(int x, int n, int oneSiteSizeInBits, unsigned int twiceS, bool isPeriodic)
{
	unsigned int mask     = (1 << oneSiteSizeInBits) - 1;
	double       sum      = 0;
	double       s        = 0.5 * twiceS;
	double       prev     = 0;
	double       val      = 0;
	double       firstVal = 0;
	for (int i = 0; i < n; ++i) {
		double tmp = (x & mask);
		assert(isValid(tmp, twiceS));
		val = tmp - s;
		if (i == 0)
			firstVal = val;
		sum += prev * val;
		prev = val;
		x >>= oneSiteSizeInBits;
	}

	if (isPeriodic)
		sum += firstVal * val;

	assert(x == 0);
	return sum;
}

int findOneSiteSizeInBits(unsigned int twiceS)
{
	for (int i = 0; i < 64; ++i) {
		unsigned int x = (1 << i);
		if (x >= twiceS)
			return i + 1;
	}

	throw std::runtime_error("N or S too big\n");
}

void f(int n, unsigned int twiceS, int targetSzPlusConst, bool isPeriodic)
{
	int    oneSiteSizeInBits = findOneSiteSizeInBits(twiceS);
	int    total             = (1 << (n * oneSiteSizeInBits));
	int    count             = 0;
	double sum               = 0;
	std::cout << "#twiceS=" << twiceS << " Total=" << total << "\n";
	std::cout << "#oneSiteSizeInBits=" << oneSiteSizeInBits << "\n";
	for (int i = 0; i < total; ++i) {
		int szPlusConst = findSzPlusConst(i, n, oneSiteSizeInBits, twiceS);
		if (szPlusConst != targetSzPlusConst)
			continue;
		double e = energy(i, n, oneSiteSizeInBits, twiceS, isPeriodic);
		sum += e;
		//	std::cout<<i<<" "<<e<<" "<<sum<<"\n";
		++count;
	}

	std::cout << (sum / count) << " " << sum << " " << count << "\n";
}

int main(int argc, char* argv[])
{
	if (argc == 1)
		return 1;

	int n          = atoi(argv[1]);
	int twiceS     = (argc >= 3) ? atof(argv[2]) : 1;
	int isPeriodic = (argc >= 4) ? atof(argv[3]) : false;
	// \sum_i (m_i + s) = Sz + s*N
	f(n, twiceS, twiceS * n / 2, isPeriodic > 0);
}
