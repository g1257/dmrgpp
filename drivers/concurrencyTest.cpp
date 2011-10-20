// BEGIN LICENSE BLOCK
/*
Copyright (c) 2011 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

------------------------------------------------------------------
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. 

Please see full open source license included in file LICENSE.
------------------------------------------------------------------

*/
// END LICENSE BLOCK

#include <string>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include "Vector.h"

typedef double FieldType;
#ifdef USE_MPI
#include "ConcurrencyMpi.h"
typedef ConcurrencyMpi<FieldType> ConcurrencyType;
#else
#include "ConcurrencySerial.h"
typedef ConcurrencySerial<FieldType> ConcurrencyType;
#endif


int main(int argc,char *argv[])
{
	ConcurrencyType concurrency(argc,argv);
	
	CommType comm1 = concurrency.newCommFromSegments(ySize);
	size_t loop1Total = atoi(argv[1]);
	concurrency.loopCreate(loop1Total,comm1);
	while(concurrency.loop(i)) {
		std::cout<<"i="<<i<<" comm1.rank="<<concurrency.rank(comm1)<<"\n";
	}
}
