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
#include <sched.h>
#include <unistd.h>
#ifdef PTHREAD_ASSIGN_AFFINITIES
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#endif
#ifdef _GNU_SOURCE
#include <errno.h>
#include <string.h>
#endif

namespace PsimagLite {

class LoadBalancerDefault {

public:

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	LoadBalancerDefault(const VectorSizeType& weights, SizeType nthreads)
	    : taskNumber_(nthreads)
	{
		VectorSizeType workLoad(nthreads, 0);

		SizeType ntasks = weights.size();
		for (SizeType i = 0; i < ntasks; ++i) {
			SizeType thread = findThreadWithLightestWork(workLoad);
			// assign work to thread
			assert(thread < taskNumber_.size());
			taskNumber_[thread].push_back(i);
			// update work loads
			assert(thread < workLoad.size());
			workLoad[thread] += weights[i];
		}

		for (SizeType i = 0; i < nthreads; ++i) {
			SizeType n = taskNumber_[i].size();
			std::cout<<n<<" Indices allocated to thread "<<i<<": ";
			for (SizeType j = 0; j < n; ++j)
				std::cout<<taskNumber_[i][j]<<" ";
			std::cout<<"\n";
		}
	}

	SizeType blockSize(SizeType threadNum) const
	{
		assert(threadNum < taskNumber_.size());
		return taskNumber_[threadNum].size();
	}

	int taskNumber(SizeType threadNum, SizeType p) const
	{
		assert(threadNum < taskNumber_.size());
		assert(p < taskNumber_[threadNum].size());
		return taskNumber_[threadNum][p];
	}

private:

	SizeType findThreadWithLightestWork(const VectorSizeType& workLoad) const
	{
		return std::min_element(workLoad.begin(), workLoad.end()) - workLoad.begin();
	}

	PsimagLite::Vector<VectorSizeType>::Type taskNumber_;
};

template<typename PthreadFunctionHolderType, typename LoadBalancerType=LoadBalancerDefault>
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

template<typename PthreadFunctionHolderType>
void *thread_function_wrapper(void *dummyPtr)
{
	PthreadFunctionStruct<PthreadFunctionHolderType> *pfs =
	        static_cast<PthreadFunctionStruct<PthreadFunctionHolderType> *>(dummyPtr);

	PthreadFunctionHolderType *pfh = pfs->pfh;

	int s = 1;
#ifdef __linux__
	s = sched_getcpu();
	if (s >= 0) pfs->cpu = s;
#endif

	SizeType blockSize = pfs->loadBalancer->blockSize(pfs->threadNum);

	for (SizeType p=0; p < blockSize; ++p) {
		SizeType taskNumber = pfs->loadBalancer->taskNumber(pfs->threadNum, p);
		if (taskNumber > pfs->total) break;
		pfh->doTask(taskNumber, pfs->threadNum);
	}

	return 0;
}

template<typename PthreadFunctionHolderType, typename LoadBalancerType=LoadBalancerDefault>
class PthreadsNg  {

public:

	typedef LoadBalancerDefault::VectorSizeType VectorSizeType;

	PthreadsNg(SizeType nPthreadsNg,int = 0)
	    : nthreads_(nPthreadsNg),cores_(1)
	{
		int cores = sysconf(_SC_NPROCESSORS_ONLN);
		cores_ = (cores > 0) ? cores : 1;
	}

	// no weights, no balancer ==> create weights, set all weigths to 1, delegate
	void loopCreate(PthreadFunctionHolderType& pfh)
	{
		SizeType ntasks = pfh.tasks();
		VectorSizeType weights(ntasks,1);
		loopCreate(pfh,weights);
	}

	// weights, no balancer ==> create balancer with weights ==> delegate
	void loopCreate(PthreadFunctionHolderType& pfh, const VectorSizeType& weights)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights, nthreads_);
		loopCreate(pfh,*loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// balancer (includes weights)
	void loopCreate(PthreadFunctionHolderType& pfh,
	                const LoadBalancerType& loadBalancer)
	{
		PthreadFunctionStruct<PthreadFunctionHolderType>* pfs;
		pfs = new PthreadFunctionStruct<PthreadFunctionHolderType>[nthreads_];
		pthread_t* thread_id = new pthread_t[nthreads_];
		pthread_attr_t** attr = new pthread_attr_t*[nthreads_];
		SizeType ntasks = pfh.tasks();

		for (SizeType j=0; j <nthreads_; j++) {
			pfs[j].pfh = &pfh;
			pfs[j].loadBalancer = &loadBalancer;
			pfs[j].threadNum = j;
			pfs[j].total = ntasks;
			pfs[j].nthreads = nthreads_;

			attr[j] = new pthread_attr_t;
			int ret = pthread_attr_init(attr[j]);
			checkForError(ret);

			setAffinity(attr[j],j,cores_);

			ret = pthread_create(&thread_id[j],
			                     attr[j],
			                     thread_function_wrapper<PthreadFunctionHolderType>,
			                     &pfs[j]);
			checkForError(ret);
		}

		for (SizeType j=0; j <nthreads_; ++j) pthread_join(thread_id[j], 0);
		for (SizeType j=0; j <nthreads_; ++j) {
			int ret = pthread_attr_destroy(attr[j]);
			checkForError(ret);
			delete attr[j];
			attr[j] = 0;
		}

		delete [] attr;

#ifndef NDEBUG
		for (SizeType j=0; j <nthreads_; j++) {
			std::cout<<"PthreadsNg: Pthread number "<<j<<" runs on core number ";
			std::cout<<pfs[j].cpu<<"\n";
		}
#endif

		delete [] thread_id;
		delete [] pfs;
	}

	String name() const { return "PthreadsNg"; }

	SizeType threads() const { return nthreads_; }

	SizeType mpiProcs() const { return 1; }

private:

	void setAffinity(pthread_attr_t* attr,
	                 SizeType threadNum,
	                 SizeType cores) const
	{
#ifdef PTHREAD_ASSIGN_AFFINITIES
		cpu_set_t* cpuset = new cpu_set_t;
		int cpu = threadNum % cores;
		CPU_ZERO(cpuset);
		CPU_SET(cpu,cpuset);
		std::size_t cpusetsize = sizeof(cpu_set_t);
		int ret = pthread_attr_setaffinity_np(attr,cpusetsize,cpuset);
		checkForError(ret);
		// clean up
		delete cpuset;
		cpuset = 0;
#endif
	}

	void checkForError(int ret) const
	{
		if (ret == 0) return;
#ifdef _GNU_SOURCE
		std::cerr<<"PthreadsNg ERROR: "<<strerror(ret)<<"\n";
#endif
	}

	SizeType nthreads_;
	SizeType cores_;
}; // PthreadsNg class

} // namespace Dmrg

/*@}*/
#endif
