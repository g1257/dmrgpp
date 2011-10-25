#include "Range.h"
#include <iostream>

typedef double RealType;

#ifndef USE_MPI
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<> ConcurrencyType;
#else
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi <RealType> ConcurrencyType;
#endif

int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	size_t total = 10;
	PsimagLite::Range < ConcurrencyType > range(0,total,concurrency);

	RealType sum = 0.0;
	while (!range.end()) {
		sum += range.index();
		range.next();
	}
	std::cout<<"sum="<<sum<<"\n";

	PsimagLite::Range < ConcurrencyType > range2(0,total,concurrency);
	for (;!range2.end();range2.next()) {
		std::cout<<"rank="<<concurrency.rank()<<" "<<range2.index()<<"\n";
	}
}
