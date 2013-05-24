/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 2.0.0]
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
/** \ingroup DMRG */
/*@{*/

/*! \file Checkpoint.h
 *
 *  checkpointing functions
 *  this class also owns the stacks since they 
 *  are so related to checkpointing
 */
#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "Stack.h"
#include "DiskStack.h"
#include "ProgressIndicator.h"
#include "ProgramGlobals.h"

namespace Dmrg {

	template<typename ParametersType,typename TargettingType>
	class Checkpoint {
	public:
		typedef typename TargettingType::RealType  RealType;
		typedef typename TargettingType::BasisWithOperatorsType BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::OperatorsType OperatorsType;
		typedef typename TargettingType::IoType IoType;
		typedef typename PsimagLite::Stack<BasisWithOperatorsType>::Type MemoryStackType;
		typedef DiskStack<BasisWithOperatorsType>  DiskStackType;

		const PsimagLite::String SYSTEM_STACK_STRING;
		const PsimagLite::String ENVIRON_STACK_STRING;

		Checkpoint(const ParametersType& parameters,SizeType rank = 0,bool debug=false) :
			SYSTEM_STACK_STRING("SystemStack"),
			ENVIRON_STACK_STRING("EnvironStack"),
			parameters_(parameters),
			enabled_(parameters_.options.find("checkpoint")!=PsimagLite::String::npos || parameters_.options.find("restart")!=PsimagLite::String::npos),
			systemDisk_(SYSTEM_STACK_STRING+parameters_.checkpoint.filename , SYSTEM_STACK_STRING+parameters_.filename,enabled_,rank),
			envDisk_(ENVIRON_STACK_STRING+parameters_.checkpoint.filename , ENVIRON_STACK_STRING+parameters_.filename,enabled_,rank),
			progress_("Checkpoint",rank)
		{
			if (!enabled_) return;
			if (parameters_.checkpoint.filename == parameters_.filename) {
				throw PsimagLite::RuntimeError("Checkpoint::ctor(...): "
						"this run will overwrite previous, throwing\n");
			}
			loadStacksDiskToMemory();
		}

		~Checkpoint()
		{
			 loadStacksMemoryToDisk();
		}

		// Not related to stacks
		void save(const BasisWithOperatorsType &pS,const BasisWithOperatorsType &pE,typename IoType::Out& io) const
		{
			PsimagLite::OstringStream msg;
			msg<<"Saving pS and pE...";
			progress_.printline(msg,std::cout);
			pS.save(io,"#CHKPOINTSYSTEM");
			pE.save(io,"#CHKPOINTENVIRON");
		}

		// Not related to stacks
		void load(BasisWithOperatorsType &pS,BasisWithOperatorsType &pE,TargettingType& psi)
		{

			typename IoType::In ioTmp(parameters_.checkpoint.filename);
			SizeType loop = ioTmp.count("#NAME=#CHKPOINTSYSTEM");
			if (loop<1) {
				std::cerr<<"There are no resumable loops in file "<<parameters_.checkpoint.filename<<"\n";
				throw PsimagLite::RuntimeError("Checkpoint::load(...)\n");
			}
			loop--;
			BasisWithOperatorsType pS1(ioTmp,"#CHKPOINTSYSTEM",loop);
			pS=pS1;
			BasisWithOperatorsType pE1(ioTmp,"#CHKPOINTENVIRON");
			pE=pE1;
			psi.load(parameters_.checkpoint.filename);
		}

		void push(const BasisWithOperatorsType &pS,const BasisWithOperatorsType &pE)
		{
			systemStack_.push(pS);
			envStack_.push(pE);
		}

		void push(const BasisWithOperatorsType &pSorE,SizeType what)
		{
			if (what==ProgramGlobals::ENVIRON) envStack_.push(pSorE);
			else systemStack_.push(pSorE);
		}

		BasisWithOperatorsType shrink(SizeType what,const TargettingType& target)
		{
			if (what==ProgramGlobals::ENVIRON) return shrink(envStack_,target);
			else return shrink(systemStack_,target);
		}

		bool operator()() const { return enabled_; }

		SizeType stackSize(SizeType what) const
		{
			if (what==ProgramGlobals::ENVIRON) return envStack_.size();
			return systemStack_.size();
		}

	private:

		//! shrink  (we don't really shrink, we just undo the growth)
		BasisWithOperatorsType shrink(MemoryStackType& thisStack,const TargettingType& target)
		{
			thisStack.pop();
			BasisWithOperatorsType& basisWithOps =  thisStack.top();
			// only updates the extreme sites:
			target.updateOnSiteForTimeDep(basisWithOps);
			return basisWithOps;
		}

		void loadStacksDiskToMemory()
		{
			PsimagLite::OstringStream msg;
			msg<<"Loading sys. and env. stacks from disk...";
			progress_.printline(msg,std::cout);

			loadStack(systemStack_,systemDisk_);
			loadStack(envStack_,envDisk_);
		}

		void loadStacksMemoryToDisk()
		{
			PsimagLite::OstringStream msg;
			msg<<"Writing sys. and env. stacks to disk...";
			progress_.printline(msg,std::cout);
			loadStack(systemDisk_,systemStack_);
			loadStack(envDisk_,envStack_);
		}

		template<typename StackType1,typename StackType2>
		void loadStack(StackType1& stackInMemory,StackType2& stackInDisk)
		{
			while (stackInDisk.size()>0) {
				BasisWithOperatorsType b = stackInDisk.top();
				stackInMemory.push(b);
				stackInDisk.pop();
			}

		}

		//! Move elsewhere
		//! returns s1+s2 if s2 has no '/',
		//! if s2 = s2a + '/' + s2b return s2a + '/' + s1 + s2b
		PsimagLite::String appendWithDir(const PsimagLite::String& s1,const PsimagLite::String& s2) const
		{
			SizeType x = s2.find("/");
			if (x==PsimagLite::String::npos) return s1 + s2;
			PsimagLite::String suf = s2.substr(x+1,s2.length());
			PsimagLite::String dir = s2.substr(0,s2.length()-suf.length());
			return dir + s1 + suf;
		}

		const ParametersType& parameters_;
		bool enabled_;
		MemoryStackType systemStack_,envStack_; // <--we're the owner
		DiskStackType systemDisk_,envDisk_;
		PsimagLite::ProgressIndicator progress_;
	}; // class Checkpoint
} // namespace Dmrg 

/*@}*/
#endif
