// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009 , UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.0.0]
[by G.A., Oak Ridge National Laboratory]

UT Battelle Open Source Software License 11242008

OPEN SOURCE LICENSE

Subject to the conditions of this License, each
contributor to this software hereby grants, free of
charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), a
perpetual, worldwide, non-exclusive, no-charge,
royalty-free, irrevocable copyright license to use, copy,
modify, merge, publish, distribute, and/or sublicense
copies of the Software.

1. Redistributions of Software must retain the above
copyright and license notices, this list of conditions,
and the following disclaimer.  Changes or modifications
to, or derivative works of, the Software should be noted
with comments and the contributor and organization's
name.

2. Neither the names of UT-Battelle, LLC or the
Department of Energy nor the names of the Software
contributors may be used to endorse or promote products
derived from this software without specific prior written
permission of UT-Battelle.

3. The software and the end-user documentation included
with the redistribution, with or without modification,
must include the following acknowledgment:

"This product includes software produced by UT-Battelle,
LLC under Contract No. DE-AC05-00OR22725  with the
Department of Energy."
 
*********************************************************
DISCLAIMER

THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
COPYRIGHT OWNER, CONTRIBUTORS, UNITED STATES GOVERNMENT,
OR THE UNITED STATES DEPARTMENT OF ENERGY BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED
STATES DEPARTMENT OF ENERGY, NOR THE COPYRIGHT OWNER, NOR
ANY OF THEIR EMPLOYEES, REPRESENTS THAT THE USE OF ANY
INFORMATION, DATA, APPARATUS, PRODUCT, OR PROCESS
DISCLOSED WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.

*********************************************************


*/
// END LICENSE BLOCK
/** \ingroup PsimagLite */
/*@{*/

/*! \file Rusage.h
 *
 * getrusage system call
 */
  
#ifndef RUSAGE_H_H
#define RUSAGE_H_H
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>

namespace PsimagLite {
	class Rusage {
	public:
		enum {USER_TIME,SYSTEM_TIME};
		typedef std::pair<size_t,size_t> PairType;

		//  who can be either RUSAGE_SELF or RUSAGE_CHILDREN
		Rusage(int who = RUSAGE_SELF) : who_(who)
		{
			getrusage(who,&rusage_);
		}
		
		void update()
		{
			getrusage(who_,&rusage_);
		}

		/* user time used or system time used */
		std::pair<size_t,size_t> time(size_t userOrSystem = USER_TIME) const
		{
			if (userOrSystem == USER_TIME)
				return PairType(
			                        rusage_.ru_utime.tv_sec,
			                        rusage_.ru_utime.tv_usec);
			return PairType(
			                rusage_.ru_stime.tv_sec,
			                rusage_.ru_stime.tv_usec);
		}

		/* maximum resident set size */
		long maxSet() const { return rusage_.ru_maxrss; }
		
		/* integral shared memory size */
		long sharedMemory() const { return rusage_.ru_ixrss; }
		
		/* integral unshared data size */
		long unsharedData() const { return rusage_.ru_idrss; }

		/* integral unshared stack size */
		long unsharedStack() const { return rusage_.ru_isrss; }

		long memory()
		{
			update();
			return sharedMemory() + unsharedData() + unsharedStack();
		}

		/* page reclaims */
		long pageReclaims() const { return rusage_.ru_minflt; }

		/* page faults */
		long pageFaults() const { return rusage_.ru_majflt; }

		/* swaps */
		long swaps() const { return rusage_.ru_nswap; }

		/* block input operations */
		long blockInputOps() const { return rusage_.ru_inblock; }

		/* block output operations */
		long blockOutputOps() const { return rusage_.ru_oublock; }

		/* messages sent */
		long messagesSent() const { return rusage_.ru_msgsnd; }

		/* messages received */
		long messagesReceived() const { return rusage_.ru_msgrcv; }

		/* signals received */
		long signalsReceived() const { return rusage_.ru_nsignals; }

		/* voluntary context switches */
		long voluntaryContextSwitches() const { return rusage_.ru_nvcsw; }
		
		/* involuntary context switches */
		long involuntaryContextSwitches() const { return rusage_.ru_nivcsw; }

	private:
		int who_;
		struct rusage rusage_;
	}; // class Rusage
} // namespace PsimagLite 

/*@}*/	
#endif // RUSAGE_H_H
