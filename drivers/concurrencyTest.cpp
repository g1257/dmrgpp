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

typedef double RealType;
#ifdef USE_MPI
#include "ConcurrencyMpi.h"
typedef PsimagLite::ConcurrencyMpi<RealType> ConcurrencyType;
#else
#include "ConcurrencySerial.h"
typedef PsimagLite::ConcurrencySerial<RealType> ConcurrencyType;
#endif


int main(int argc,char *argv[])
{
	if (argc<3) {
		std::string s(argv[0]);
		s += ": Needs loopTotal and segmentSize as args\n";
		throw std::runtime_error(s.c_str());
	}
	typedef ConcurrencyType::CommType CommType;
	ConcurrencyType concurrency(argc,argv);
	
	
	size_t loop1Total = atoi(argv[1]);
	size_t ySize = atoi(argv[2]);
	CommType comm1 = concurrency.newCommFromSegments(ySize);
	
	concurrency.loopCreate(loop1Total,comm1);
	size_t i = 0;
	while(concurrency.loop(i)) {
		std::cout<<"i="<<i<<" comm1.rank="<<concurrency.rank(comm1)<<"\n";
	}
}
