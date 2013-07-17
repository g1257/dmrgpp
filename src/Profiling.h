/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

// A class to profile a scope of code
#ifndef PROFILING_H_
#define PROFILING_H_

#include <iostream>
#include "MemoryUsage.h"
#include "String.h"

namespace PsimagLite {

std::ostream& operator<<(std::ostream& os,const std::pair<SizeType,SizeType>& p)
{

	os<<p.first<<" "<<p.second<<" ";
	return os;
}

class  Profiling {

	double diff(double x1,double x2) const
	{
		return x1 - x2;
	}

public:

	Profiling(const String& s,std::ostream& os = std::cout)
	    : message_(s),
	      memoryUsage_("/proc/self/stat"),
	      start_(memoryUsage_.time()),
	      isDead_(false),
	      os_(os)
	{
		os_<<"Profiling: Starting clock for "<<s<<"\n";
	}

	~Profiling()
	{
		killIt();
	}

	void end()
	{
		killIt();
	}

private:

	void killIt()
	{
		if (isDead_) return;
		double end = memoryUsage_.time();
		double elapsed = diff(end,start_);
		os_<<"Profiling: Stoping clock for "<<message_;
		os_<<" elapsed="<<elapsed<<"\n";
		std::cout<<" start="<<start_<<" end="<<end<<"\n";
		isDead_ = true;
	}

	String message_;
	MemoryUsage memoryUsage_;
	double start_;
	bool isDead_;
	std::ostream& os_;
}; // Profiling
} // PsimagLite

#endif

