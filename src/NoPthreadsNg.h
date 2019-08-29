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

/*! \file NoPthreadsNg .h
 *
 *  A C++ NoPthreadsNg class that implements the Concurrency interface
 *
 */
#ifndef PSI_NOPTHREADS_NG_H
#define PSI_NOPTHREADS_NG_H
#include "LoadBalancerDefault.h"
#include "CodeSectionParams.h"

namespace PsimagLite {

template<typename PthreadFunctionHolderType, typename LoadBalancerType=LoadBalancerDefault>
class NoPthreadsNg  {

public:

	typedef LoadBalancerDefault::VectorSizeType VectorSizeType;

	NoPthreadsNg(const CodeSectionParams& cs)
	{
		if (cs.npthreads != 1)
			throw PsimagLite::RuntimeError("NoPthreadsNg: ctor with threads != 1\n");
	}

	bool affinities() const { return false; }

	size_t stackSize() const { return 0; }

	// no weights, no balancer ==> create weights, set all weigths to 1, delegate
	void loopCreate(PthreadFunctionHolderType& pfh)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(pfh.tasks(), 1);
		loopCreate(pfh, *loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// weights, no balancer ==> create balancer with weights ==> delegate
	void loopCreate(PthreadFunctionHolderType& pfh, const VectorSizeType& weights)
	{
		LoadBalancerType* loadBalancer = new LoadBalancerType(weights, 1);
		loopCreate(pfh,*loadBalancer);
		delete loadBalancer;
		loadBalancer = 0;
	}

	// balancer (includes weights)
	void loopCreate(PthreadFunctionHolderType& pfh,
	                const LoadBalancerType&)
	{
		SizeType ntasks = pfh.tasks();
		for (SizeType taskNumber = 0; taskNumber < ntasks; ++taskNumber) {
			pfh.doTask(taskNumber, 0);
		}
	}

	String name() const { return "NoPthreadsNg"; }

	SizeType threads() const { return 1; }

	SizeType mpiProcs() const { return 1; }

	void setAffinities(bool) {}
}; // NoPthreadsNg class

} // namespace Dmrg

/*@}*/
#endif
