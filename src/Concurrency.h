/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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
/** \ingroup PsimagLite */
/*@{*/

/*! \file Concurrency.h
 *
 */
#ifndef CONCURRENCY_HEADER_H
#define CONCURRENCY_HEADER_H
#include <stdexcept>
#include <cassert>
#include <sys/syscall.h>
#include <unistd.h>
#include "Vector.h"
#include "Mpi.h"
#include "FloatingPoint.h"
#include "LAPACK.h"

namespace PsimagLite {

class Concurrency {

	typedef MpiDisabled MpiDisabledType;

public:

	static SizeType mode;
	static SizeType npthreads;
	static bool setAffinitiesDefault;

#ifndef USE_PTHREADS
	typedef int MutexType;

	static void mutexLock(MutexType*)
	{}

	static void mutexUnlock(MutexType*)
	{}

	static void mutexInit(MutexType*)
	{}

	static void mutexDestroy(MutexType*)
	{}

#else

	#include <pthread.h>
	typedef pthread_mutex_t MutexType;

	static void mutexInit(MutexType* mutex)
	{
		pthread_mutex_init(mutex, 0);
	}

	static void mutexDestroy(MutexType* mutex)
	{
		pthread_mutex_destroy(mutex);
	}

	static void mutexLock(MutexType* mutex)
	{
		pthread_mutex_lock(mutex);
	}

	static void mutexUnlock(MutexType* mutex)
	{
		pthread_mutex_unlock(mutex);
	}

#endif

#ifndef USE_MPI
	typedef int CommType;
#else
	typedef MPI::CommType CommType;
#endif

	enum {SERIAL=0,PTHREADS=1,MPI=2,PTHREADS_AND_MPI=3};

	static SizeType storageSize(SizeType npthreads)
	{
		switch (mode) {
		case SERIAL:
			assert(npthreads == 1);
		case PTHREADS:
		case PTHREADS_AND_MPI:
			return npthreads;
		case MPI:
			return 1;
		}
		throw RuntimeError("storageSize: wrong mode\n");
	}

	static SizeType storageIndex(SizeType threadNum)
	{
		switch (mode) {
		case SERIAL:
			assert(threadNum == 0);
		case PTHREADS:
		case PTHREADS_AND_MPI:
			return threadNum;
		case MPI:
			return 0;
		}
		throw RuntimeError("storageIndex: wrong mode\n");
	}

#ifndef USE_MPI
	Concurrency(int*, char ***,size_t nthreads)
#else
	Concurrency(int* argc, char *** argv,size_t nthreads)
#endif
	{
		FloatingPoint::enableExcept();
		npthreads = nthreads;
		mode = 0;
#ifdef USE_PTHREADS
		mode |= 1;
		if (!psimag::LAPACK::isThreadSafe())
			std::cerr<<"WARNING: You LAPACK might not be thread safe\n";
#else
		if (nthreads != 1)
			throw RuntimeError("nthreads>1 but no USE_PTHREADS support compiled\n");
#endif
#ifdef USE_MPI
		MPI::init(argc,argv);
		mode |= 2;
#endif
	}

	~Concurrency()
	{
		MPI::finalize();
	}

	static bool root(MPI::CommType comm = MPI::COMM_WORLD)
	{
		return (MPI::commRank(comm) == 0);
	}

	static SizeType nprocs(MPI::CommType comm = MPI::COMM_WORLD)
	{
		return MPI::commSize(comm);
	}

	static SizeType rank(MPI::CommType comm = MPI::COMM_WORLD)
	{
		return MPI::commRank(comm);
	}

	static bool hasMpi()
	{
		return (mode & MPI);
	}

	static bool hasPthreads()
	{
		return (mode & PTHREADS);
	}

	static void mpiDisable(String label)
	{
		if (!hasMpi()) return;
		mpiDisabled_.disable(label);
	}

	static void mpiDisableIfNeeded(SizeType& mpiRank,
	                               SizeType& blockSize,
	                               String label,
	                               SizeType total)
	{
		if (!hasMpi()) return;
		if (!mpiDisabled_(label)) return;
		mpiRank = 0;
		blockSize = total;
		if (!hasPthreads()) return;
		String str(__FILE__);
		str += " mpiDisableIfNeeded label = " + label + "\n";
		throw RuntimeError(str);
	}

	static bool isMpiDisabled(String label)
	{
		if (!hasMpi()) return false;
		return mpiDisabled_(label);
	}

	static void setOptions(SizeType nthreads, bool affinities)
	{
		npthreads = nthreads;
		setAffinitiesDefault = affinities;
		if (nthreads==1) return;
#ifndef USE_PTHREADS
		PsimagLite::String message1(__FILE__);
		message1 += " FATAL: You are requesting nthreads>0 but you ";
		message1 += "did not compile with USE_PTHREADS enabled\n";
		message1 += " Either set Threads=1 in the input file (you won't ";
		message1 += "have threads though) or\n";
		message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile ";
		message1 += "and recompile\n";
		throw PsimagLite::RuntimeError(message1.c_str());
#else
		std::cout<<"Concurrency::npthreads="<<npthreads<<"\n";
		std::cout<<"Concurrency::setAffinitiesDefault="<<setAffinitiesDefault<<"\n";
#endif
	}

private:

	static MpiDisabledType mpiDisabled_;
};

} // namespace PsimagLite

/*@}*/
#endif

