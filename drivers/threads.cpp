/* Copyright (c) 2009-2013, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
#include "ConcurrencySerial.h"
#include <iostream>
#include <cstdlib>
#ifdef USE_PTHREADS
#define PTHREADS_NAME PsimagLite::Pthreads
#include "Pthreads.h"
#else
#define PTHREADS_NAME PsimagLite::NoPthreads
#include "NoPthreads.h"
#endif


class MyHelper {

public:

	typedef double RealType;

	void thread_function_(size_t threadNum,size_t blockSize,size_t total,pthread_mutex_t* myMutex)
	{
		for (size_t p=0;p<blockSize;p++) {
			size_t taskNumber = threadNum*blockSize + p;
			if (taskNumber>=total) break;
			std::cout<<"This is thread number "<<threadNum<<" and taskNumber="<<taskNumber<<"\n";
		}
	}

}; // class MyHelper

typedef PsimagLite::ConcurrencySerial<double> ConcurrencyType;
typedef MyHelper HelperType;

int main(int argc,char *argv[])
{
	if (argc!=3) {
		std::cout<<"USAGE: "<<argv[0]<<" nthreads ntasks\n";
		return 1;
	}

	size_t nthreads  = atoi(argv[1]);
	size_t ntasks = atoi(argv[2]);

	std::cout<<"Using "<<nthreads<<" threads.\n";

	ConcurrencyType concurrency(argc,argv);

	PTHREADS_NAME<HelperType> threadObject;

	PTHREADS_NAME<HelperType>::setThreads(nthreads);

	HelperType helper;

	threadObject.loopCreate(ntasks,helper,concurrency);

}

