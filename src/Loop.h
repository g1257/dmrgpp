// BEGIN LICENSE BLOCK
/*
Copyright (c) 2011, UT-Battelle, LLC
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

/*! \file Loop.h
 *
 * A loop class that can be parallelized
 */
#ifndef LOOP_HEADER_H
#define LOOP_HEADER_H
#include <stdexcept>
#include <vector>
#include "ConcurrencySerial.h" // for the default Loop

namespace PsimagLite {

	template<typename ConcurrencyType>
	class Loop {

	public:

		typedef typename ConcurrencyType::CommType CommType;
		static const CommType COMM_WORLD;

		Loop(ConcurrencyType& concurrency,
		     size_t total,
		     const std::vector<size_t>& weights,
		     CommType mpiComm=COMM_WORLD)
		: concurrency_(concurrency),
		  total_(total),
		  step_(0),
		  nprocs_(concurrency.nprocs(mpiComm)),
		  indicesOfThisProc_(nprocs_),
		  assigned_(false)
		{
			init(weights,mpiComm);
		}

		Loop(ConcurrencyType& concurrency,
		                size_t total,
		                CommType mpiComm=COMM_WORLD)
		: concurrency_(concurrency),
		  total_(total),
		  step_(0),
		  nprocs_(concurrency.nprocs(mpiComm)),
		  indicesOfThisProc_(nprocs_),
		  assigned_(false)
		{
			std::vector<size_t> weights(total,1);
			init(weights,mpiComm);
		}

		bool next(size_t &i)
		{
			if (!assigned_) return false;
			
			if (step_>=myIndices_.size())  {
				step_ = 0;
				return false;
			}

			i=myIndices_[step_];

			if (i>=total_ ) {
				step_= 0; 
				return false;
			}
			step_++;
			return true;
		}

	private:

		ConcurrencyType& concurrency_;
		size_t total_; // total number of indices total_(total),
		size_t step_; // step within this processor
		size_t nprocs_;
		std::vector<std::vector<size_t> > indicesOfThisProc_; // given rank and step it maps the index
		bool assigned_;
		std::vector<size_t> myIndices_; // indices assigned to this processor

		void init(const std::vector<size_t>& weights,CommType mpiComm)
		{
			size_t r1=concurrency_.rank(mpiComm);
			
			// distribute the load among the processors
			std::vector<size_t> loads(nprocs_,0);
			
			for (size_t i=0;i<total_;i++) {
				size_t r = findLowestLoad(loads);
				indicesOfThisProc_[r].push_back(i);
				loads[r] += weights[i];
				if (r==r1) assigned_=true;
			}
			// set myIndices_
			myIndices_=indicesOfThisProc_[r1];
			//MPI_Barrier(mpiComm);
		}

		size_t findLowestLoad(std::vector<size_t> const &loads) const
		{
			size_t x= 1000000;
			size_t ret=0;
			for (size_t i=0;i<loads.size();i++) {
				if (loads[i]<x) {
					x=loads[i];
					ret =i;
				}
			}
			return ret;
		}
		
	}; // class Loop

	template<typename ConcurrencyType>
	const typename Loop<ConcurrencyType>::CommType 
	Loop<ConcurrencyType>::COMM_WORLD = ConcurrencyType::COMM_WORLD;


	//! Specialization for performance reasons
	template<>
	class Loop<ConcurrencySerial<> > {
	
		typedef ConcurrencySerial<> ConcurrencyType;

	public:

		typedef typename ConcurrencyType::CommType CommType;
		Loop(ConcurrencyType& concurrency,
		     size_t total,
		     std::vector<size_t> const &weights,
		     CommType mpiComm=0)
		: total_(total),step_(0)
		{}
		
		Loop(ConcurrencyType& concurrency,
		     size_t total,
		     CommType mpiComm=0)
		: total_(total),step_(0)
		{}

		bool next(size_t &i)
		{
			i=step_;
			
			if (i>=total_ ) {
				step_= 0; 
				return false;
			}
			step_++;
			return true;
		}

	private:
		size_t total_; // total number of indices total_(total),
		size_t step_; // step within this processor
	}; // class Loop
} // namespace Dmrg

/*@}*/	
#endif // LOOP_HEADER_H
