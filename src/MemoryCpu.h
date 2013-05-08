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
/** \ingroup PsimagLite */
/*@{*/

/*!
 *
 *
 */
#ifndef MEMORY_CPU_H
#define MEMORY_CPU_H

#include <vector>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include "ProgressIndicator.h"

class MemoryCpu {

public:

	MemoryCpu()
	    : progress_("MemoryCpu",0)
	{}

	void deallocate(void *p)
	{
		assert(p);
		free(p);
		p=0;
		std::ostringstream msg;
		msg<<"Freed "<<p;
		progress_.printline(msg,std::cout);
	}

	void* allocate(size_t x)
	{
		void *p = malloc(x);
		std::ostringstream msg;
		msg<<"Allocated starting at "<<p<<" "<<x<<" bytes";
		progress_.printline(msg,std::cout);
		return p;
	}

private:

	PsimagLite::ProgressIndicator progress_;

}; // class MemoryCpu

MemoryCpu globalMemoryCpu;

/*@}*/
#endif // MEMORY_CPU_H
