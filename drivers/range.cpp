/** How to compile and run this driver
 * 
 * Serial version:
 * 
 * g++ -g3 -DNDEBUG  -Werror -Wall -I../src -I../src/JSON \
 * -I../src/JSON/JsonParser -lm  -lpthread   range.cpp   -o range \
 * /usr/lib64/libblas.so.3 /usr/lib64/liblapack.so.3
 * 
 * And run it with:
 * 
 * ./range
 * 
 * Parallel version:
 * 
 * mpicxx -DUSE_MPI -g3 -DNDEBUG  -Werror -Wall -I../src -I../src/JSON \
 * -I../src/JSON/JsonParser -lm  -lpthread   range.cpp   -o range \
 * /usr/lib64/libblas.so.3 /usr/lib64/liblapack.so.3
 * 
 * And run it with:
 * 
 * your batch system script
 *
 */

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
	concurrency.reduce(sum);
	if (concurrency.root()) std::cout<<"sum="<<sum<<"\n";

	PsimagLite::Range < ConcurrencyType > range2(0,total,concurrency);
	for (;!range2.end();range2.next()) {
		std::cout<<"rank="<<concurrency.rank()<<" "<<range2.index()<<"\n";
	}
}
