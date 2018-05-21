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
#include "PsimagLite.h"
#include "ProgressIndicator.h"

namespace PsimagLite {

class  Profiling {

	double diff(double x1,double x2) const
	{
		return x1 - x2;
	}

public:

	Profiling(String caller, std::ostream& os)
	    : progressIndicator_(caller),
	      memoryUsage_("/proc/self/stat"),
	      isDead_(false),
	      os_(os)
	{
		OstringStream msg;
		msg << "starting clock";
		progressIndicator_.printline(msg, os);
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
		OstringStream msg;
		msg << "stopping clock";
		progressIndicator_.printline(msg, os_);
		isDead_ = true;
	}

	ProgressIndicator progressIndicator_;
	MemoryUsage memoryUsage_;
	bool isDead_;
	std::ostream& os_;
}; // Profiling
} // PsimagLite

#endif

