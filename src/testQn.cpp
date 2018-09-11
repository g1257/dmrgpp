#include "NotReallySort.h"
#include <cstdlib>
#include "ProgramGlobals.h"
#include <tr1/unordered_map>

typedef Dmrg::Qn QnType;
typedef QnType::VectorQnType VectorQnType;
typedef QnType::VectorSizeType VectorSizeType;
typedef QnType::PairSizeType PairSizeType;

void randomSzPlusConst(VectorSizeType& szPlusConst)
{
	SizeType nup = static_cast<int>(20*drand48());
	SizeType ndown = static_cast<int>(20*drand48());
	assert(szPlusConst.size() == 2);
	szPlusConst[0] = nup + ndown;
	szPlusConst[1] = nup;
}

void randomQn(VectorQnType& qn, VectorSizeType& szPlusConst)
{
	SizeType n = qn.size();
	for (SizeType i = 0; i < n; ++i) {
		randomSzPlusConst(szPlusConst);
		bool odd = szPlusConst[0] & 1;
		qn[i] = QnType(odd, szPlusConst, PairSizeType(0,0), 0);
	}
}


int main(int argc, char **argv)
{
	if (argc != 2) {
		std::cerr<<"Needs number of qns\n";
		return 1;
	}

	SizeType n = atoi(argv[1]);
	VectorSizeType szPlusConst(2, 0);
	VectorQnType qns(n, QnType(false, szPlusConst, PairSizeType(0,0), 0));
	randomQn(qns, szPlusConst);
	VectorSizeType outNumber;
	VectorQnType outQns;
	VectorSizeType offset;
	VectorSizeType inNumbers(n, 0);
	for (SizeType i = 0; i < n; ++i) inNumbers[i] = i;

	 Dmrg::NotReallySort notReallySort;
	 notReallySort(outNumber, outQns, offset, inNumbers, qns, Dmrg::ProgramGlobals::VERBOSE_YES);
}


