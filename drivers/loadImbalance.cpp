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
#define USE_PTHREADS_OR_NOT_NG
#include "Parallelizer.h"

class MyHelper {

	typedef PsimagLite::Concurrency ConcurrencyType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

public:

	MyHelper(SizeType ntasks, SizeType nthreads)
	    : weight_(ntasks),
	      x_(nthreads,0)
	{
		srand48(1234);
		for (SizeType i = 0; i < ntasks; ++i) {
			double x = 10*drand48();
			weight_[i] = 1 + static_cast<SizeType>(x);
			std::cout<<weight_[i]<<" ";
		}

		std::cout<<"\n";
	}

	SizeType tasks() const { return weight_.size(); }

	int result() const
	{
		return x_[0];
	}

	const VectorSizeType& weights() { return weight_; }

	void doTask(SizeType taskNumber, SizeType threadNum)
	{
		for (SizeType i = 0; i < weight_[taskNumber]; ++i)
			x_[threadNum] += (taskNumber + i);
	}

	void sync()
	{
		for (SizeType i = 1; i < x_.size(); ++i)
			x_[0] += x_[i];
	}

private:

	VectorSizeType weight_;
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
	ParallelizerType threadObject(ConcurrencyType::codeSectionParams);

	HelperType helper(ntasks, nthreads);

	std::cout<<"Using "<<threadObject.name();
	std::cout<<" with "<<threadObject.threads()<<" threads.\n";
	threadObject.loopCreate(helper, helper.weights());
	helper.sync();
	std::cout<<"Sum of all tasks= "<<helper.result()<<"\n";
}

