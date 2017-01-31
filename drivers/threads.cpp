/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/
#include "Concurrency.h"
#include <iostream>
#include <cstdlib>
#include "Parallelizer.h"

class MyHelper {

	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	MyHelper(SizeType ntasks, SizeType nthreads)
	    : ntasks_(ntasks), x_(nthreads,0)
	{}

	SizeType tasks() const { return ntasks_; }

	int result() const
	{
		return x_[0];
	}

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		x_[threadNum] += taskNumber;
	}

	void sync()
	{
		for (SizeType i = 1; i < x_.size(); ++i)
			x_[0] += x_[i];
	}

private:

	SizeType ntasks_;
	VectorSizeType x_;
}; // class MyHelper


int main(int argc,char *argv[])
{
	typedef PsimagLite::Concurrency ConcurrencyType;


	if (argc!=3) {
		std::cout<<"USAGE: "<<argv[0]<<" nthreads ntasks\n";
		return 1;
	}

	SizeType nthreads  = atoi(argv[1]);
	SizeType ntasks = atoi(argv[2]);

	ConcurrencyType concurrency(&argc,&argv,nthreads);

	typedef MyHelper HelperType;
	typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);

	HelperType helper(ntasks, nthreads);

	std::cout<<"Using "<<threadObject.name();
	std::cout<<" with "<<threadObject.threads()<<" threads.\n";
	threadObject.loopCreate(helper);
	helper.sync();
	std::cout<<"Sum of all tasks= "<<helper.result()<<"\n";
}

