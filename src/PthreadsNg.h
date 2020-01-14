/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[PsimagLite, Version 1.]
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

---------------------------------------------------------------
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

---------------------------------------------------------------

*/
/** \ingroup PsimagLite */
/*@{*/

/*! \file PthreadsNg .h
 *
 *  A C++ PthreadsNg class that implements the Concurrency interface
 *
 */
#ifndef PSI_PTHREADS_NG_H
#define PSI_PTHREADS_NG_H

#include <pthread.h>
#include <iostream>
#include <algorithm>
#include "Vector.h"
#include "LoadBalancerDefault.h"
#include <sched.h>
#include <unistd.h>
#include "TypeToString.h"
#include "CodeSectionParams.h"

#ifdef __linux__
#include <sys/types.h>
#include <sys/syscall.h>
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#ifdef _GNU_SOURCE
#include <errno.h>
#include <string.h>
#endif
#ifndef PTHREAD_STACK_MIN
#define PTHREAD_STACK_MIN 16384
#endif

namespace PsimagLite {

template<typename PthreadFunctionHolderType, typename LoadBalancerType>
struct PthreadFunctionStruct {
	PthreadFunctionStruct()
	    : pfh(0),loadBalancer(0),threadNum(0),nthreads(0),total(0),cpu(0)
	{}

	PthreadFunctionHolderType* pfh;
	const LoadBalancerType *loadBalancer;
	int threadNum;
	SizeType nthreads;
	SizeType total;
	SizeType cpu;
};

template<typename PthreadFunctionHolderType, typename LoadBalancerType>
void *thread_function_wrapper(void *dummyPtr)
{
	PthreadFunctionStruct<PthreadFunctionHolderType, LoadBalancerType> *pfs =
	        static_cast<PthreadFunctionStruct<PthreadFunctionHolderType, LoadBalancerType> *>
	        (dummyPtr);

	PthreadFunctionHolderType *pfh = pfs->pfh;

	int s = 0;
#ifdef __linux__
	s = sched_getcpu();
#endif
	if (s >= 0) pfs->cpu = s;

	SizeType blockSize = pfs->loadBalancer->blockSize(pfs->threadNum);

	for (SizeType p=0; p < blockSize; ++p) {
		SizeType taskNumber = pfs->loadBalancer->taskNumber(pfs->threadNum, p);
		if (taskNumber >= pfs->total) break;
		pfh->doTask(taskNumber, pfs->threadNum);
	}

	int retval = 0;
	pthread_exit(static_cast<void*>(&retval));
	return 0;
}

template<typename PthreadFunctionHolderType, typename LoadBalancerType=LoadBalancerDefault>
class PthreadsNg  {

	static const int MAX_CPUS = 1024;

public:

	typedef LoadBalancerDefault::VectorSizeType VectorSizeType;

	PthreadsNg(const CodeSectionParams& codeSectionParams)
	    : nthreads_(codeSectionParams.npthreads),
	      setAffinities_(codeSectionParams.setAffinities),
	      stackSize_(codeSectionParams.stackSize)
	{}

	bool affinities() const { return setAffinities_; }

	size_t stackSize() const { return stackSize_; }

	void setAffinities(bool flag)
	{
		setAffinities_ = flag;
	}

	// no weights, no balancer ==> create weights, set all weigths to 1, delegate
	void loopCreate(PthreadFunctionHolderType& pfh)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(pfh.tasks(), nthreads_);
		loopCreate(pfh, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// weights, no balancer ==> create balancer with weights ==> delegate
	void loopCreate(PthreadFunctionHolderType& pfh, const VectorSizeType& weights)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights, nthreads_);
		loopCreate(pfh, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// balancer (includes weights)
	void loopCreate(PthreadFunctionHolderType& pfh,
	                const LoadBalancerType& loadBalancer)
	{
		SizeType ntasks = pfh.tasks();
		SizeType actualThreads = std::min(nthreads_, ntasks);
		PthreadFunctionStruct<PthreadFunctionHolderType, LoadBalancerType>* pfs;
		pfs = new PthreadFunctionStruct<PthreadFunctionHolderType, LoadBalancerType>[actualThreads];
		pthread_t* thread_id = new pthread_t[actualThreads];
		pthread_attr_t** attr = new pthread_attr_t*[actualThreads];

#ifndef __APPLE__
		cpu_set_t cpuset;


		if (setAffinities_) {
			static bool firstCall = true;
			pid_t pid = getpid(); // always successfull
			getPidAffinity(&cpuset, pid);
			int ret = sched_setaffinity(pid, sizeof(cpu_set_t), &cpuset);
			checkForError(ret);
			if (ret == 0 && firstCall) {
				printAffinity(pid, &cpuset);
				firstCall = false;
			}
		}
#endif

		for (SizeType j=0; j <actualThreads; j++) {
			pfs[j].pfh = &pfh;
			pfs[j].loadBalancer = &loadBalancer;
			pfs[j].threadNum = j;
			pfs[j].total = ntasks;
			pfs[j].nthreads = actualThreads;

			attr[j] = new pthread_attr_t;
			int ret = (stackSize_ > 0) ? pthread_attr_setstacksize(attr[j], stackSize_) : 0;
			if (ret != 0) {
				std::cerr<<__FILE__;
				std::cerr<<"\tpthread_attr_setstacksize() has returned non-zero "<<ret<<"\n";
				std::cerr<<"\tIt is possible (but no certain) that the following error";
				std::cerr<<"\thappened.\n";
				std::cerr<<"\tEINVAL The stack size is less than ";
				std::cerr<<"PTHREAD_STACK_MIN (16384) bytes.\n";
				std::cerr<<"\tI will ignore this error and let you continue\n";
			}

			ret = pthread_attr_init(attr[j]);
			checkForError(ret);

#ifndef __APPLE__
			if (setAffinities_)
				setAffinity(attr[j], &cpuset, j);
#endif

			ret = pthread_create(&thread_id[j],
			                     attr[j],
			                     thread_function_wrapper<PthreadFunctionHolderType,
			                     LoadBalancerType>,
			                     &pfs[j]);
			checkForError(ret);
		}

		for (SizeType j=0; j <actualThreads; ++j) pthread_join(thread_id[j], 0);
		for (SizeType j=0; j <actualThreads; ++j) {
			int ret = pthread_attr_destroy(attr[j]);
			checkForError(ret);
			delete attr[j];
			attr[j] = 0;
		}

		delete [] attr;
		delete [] thread_id;
		delete [] pfs;
	}

	String name() const { return "PthreadsNg"; }

	//SizeType threads() const { return nthreads_; }

	SizeType mpiProcs() const { return 1; }

private:

#ifndef __APPLE__
	void getPidAffinity(cpu_set_t* cpuset, pid_t pid) const
	{

		CPU_ZERO(cpuset);
		int ret = sched_getaffinity(pid, sizeof(cpu_set_t), cpuset);
		checkForError(ret);
		if (ret != 0) {
			CPU_ZERO(cpuset);
			return;
		}
	}

	void setAffinity(pthread_attr_t* attr, cpu_set_t* cpuset, SizeType threadNum) const
	{
		// pick a cpu from cpuset
		int chosenCpu = getOneCpuFromSet(cpuset);
		if (chosenCpu < 0) {
			std::cerr<<"setAffinity: no cpus left in set for thread "<<threadNum<<"\n";
			return;
		}

		cpu_set_t mynewset;
		CPU_ZERO(&mynewset);
		CPU_SET(chosenCpu, &mynewset); // add cpu to new set
		int ret = pthread_attr_setaffinity_np(attr, sizeof(cpu_set_t), &mynewset);
		checkForError(ret);
		if (ret != 0)
			return;

		// remove cpu from set
		CPU_CLR(chosenCpu, cpuset);
	}

	int getOneCpuFromSet(cpu_set_t* cpuset) const
	{
		int count = CPU_COUNT(cpuset);
		if (count == 0) return -1;

		int chosenCpu = -1;
		for (int i = 0; i < MAX_CPUS; ++i) {
			if (CPU_ISSET(i, cpuset) == 0) continue;
			chosenCpu = i;
			break;
		}

		return chosenCpu;
	}

	void printAffinity(pid_t pid, cpu_set_t* cpuset) const
	{
		std::cout<<"pid "<<pid<<" sched_setaffinity with count ";
		std::cout<<CPU_COUNT(cpuset)<<": ";
		for (int i = 0; i < MAX_CPUS; ++i) {
			if (CPU_ISSET(i, cpuset) == 0) continue;
			std::cout<<i<<" ";
		}

		std::cout<<"\n";
	}
#endif

	void checkForError(int ret) const
	{
		if (ret == 0) return;
		std::cerr<<"PthreadsNg ERROR: "<<strerror(ret)<<"\n";
	}

	SizeType nthreads_;
	bool setAffinities_;
	size_t stackSize_;
}; // PthreadsNg class

} // namespace Dmrg

/*@}*/
#endif
