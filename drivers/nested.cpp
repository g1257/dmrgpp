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
#include "../src/Concurrency.h"
#include <iostream>
#include <cstdlib>
#define USE_PTHREADS_OR_NOT_NG
#include "../src/Parallelizer.h"

class InnerHelper {

	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	InnerHelper(SizeType ntasks, SizeType nthreads, int outerTask)
	    : ntasks_(ntasks), x_(nthreads,0), offset_(outerTask*ntasks)
	{}

	SizeType tasks() const { return ntasks_; }

	int result() const
	{
		return x_[0];
	}

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		x_[threadNum] += taskNumber + offset_;
	}

	void sync()
	{
		for (SizeType i = 1; i < x_.size(); ++i)
			x_[0] += x_[i];
	}

private:

	SizeType ntasks_;
	VectorSizeType x_;
	SizeType offset_;
}; // class InnerHelper

class MyHelper {

	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	MyHelper(SizeType ntasks,
	         SizeType nthreadsOuter,
	         SizeType ntasksInner,
	         SizeType nthreadsInner)
	    : ntasks_(ntasks),
	      x_(nthreadsOuter,0),
	      ntasksInner_(ntasksInner),
	      nthreadsInner_(nthreadsInner)
	{}

	SizeType tasks() const { return ntasks_; }

	int result() const
	{
		return x_[0];
	}

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		typedef PsimagLite::Parallelizer<InnerHelper> ParallelizerType;
		ParallelizerType threadObject(nthreadsInner_,
		                              PsimagLite::MPI::COMM_WORLD);
		InnerHelper helper(ntasksInner_, nthreadsInner_, taskNumber);
		threadObject.loopCreate(helper);
		helper.sync();
		x_[threadNum] += helper.result();
	}

	void sync()
	{
		for (SizeType i = 1; i < x_.size(); ++i)
			x_[0] += x_[i];
	}

private:

	SizeType ntasks_;
	VectorSizeType x_;
	SizeType ntasksInner_;
	SizeType nthreadsInner_;
}; // class MyHelper


int main(int argc,char *argv[])
{
	typedef PsimagLite::Concurrency ConcurrencyType;


	if (argc != 5) {
		std::cout<<"USAGE: "<<argv[0]<<" threadsOuter tasksOuter threadsInner tasksInner\n";
		return 1;
	}

	SizeType nthreadsOuter  = atoi(argv[1]);
	SizeType ntasks = atoi(argv[2]);
	SizeType nthreadsInner  = atoi(argv[3]);
	SizeType ntasksInner = atoi(argv[4]);

	ConcurrencyType concurrency(&argc,&argv,1);

	typedef MyHelper HelperType;
	typedef PsimagLite::Parallelizer<HelperType> ParallelizerType;
	ParallelizerType threadObject(nthreadsOuter,
	                              PsimagLite::MPI::COMM_WORLD);

	HelperType helper(ntasks, nthreadsOuter, ntasksInner, nthreadsInner);

	std::cout<<"Using "<<threadObject.name();
	std::cout<<" with "<<threadObject.threads()<<" threads.\n";
	threadObject.loopCreate(helper);
	helper.sync();
	std::cout<<"Sum of all tasks= "<<helper.result()<<"\n";
}

