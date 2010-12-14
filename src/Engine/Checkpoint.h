// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file Checkpoint.h
 *
 *  checkpointing functions
 */
#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include "Utils.h"

namespace Dmrg {

	template<typename RealType,typename BasisType,typename ParametersType,typename StackType,typename IoType> 
	class Checkpoint {
		public:
			enum {SYSTEM,ENVIRON};
			
			const std::string SYSTEM_STACK_STRING;
			const std::string ENVIRON_STACK_STRING;
			
			Checkpoint(const ParametersType& parameters,size_t rank = 0,bool debug=false) :
				SYSTEM_STACK_STRING("SystemStack"),
				ENVIRON_STACK_STRING("EnvironStack"),
				parameters_(parameters),
				enabled_(parameters_.options.find("checkpoint")!=std::string::npos),
				systemStack_((enabled_) ? SYSTEM_STACK_STRING+parameters_.checkpoint.filename : SYSTEM_STACK_STRING+parameters_.filename,enabled_,rank),
				envStack_((enabled_) ? ENVIRON_STACK_STRING+parameters_.checkpoint.filename : ENVIRON_STACK_STRING+parameters_.filename,enabled_,rank)
			{
			}
			
			void save(const BasisType &pS,const BasisType &pE,size_t loop,typename IoType::Out& io) const
			{
				pS.save(io,"#CHKPOINTSYSTEM");
				pE.save(io,"#CHKPOINTENVIRON");
			}

			void load(BasisType &pS,BasisType &pE)
			{
				typename IoType::In ioTmp(parameters_.checkpoint.filename);
				BasisType pS1(ioTmp,"#CHKPOINTSYSTEM",parameters_.checkpoint.index);
				pS=pS1;
				BasisType pE1(ioTmp,"#CHKPOINTENVIRON");
				pE=pE1;
			}

			void push(const BasisType &pS,const BasisType &pE)
			{	
				systemStack_.push(pS);
				envStack_.push(pE);
			}

			void push(const BasisType &pSorE,size_t what)
			{
				if (what==ENVIRON) envStack_.push(pSorE);
				else systemStack_.push(pSorE);
			}
			
			void shrink(BasisType &pSorE,size_t what)
			{
				if (what==ENVIRON) shrink(pSorE,envStack_);
				else shrink(pSorE,systemStack_);
			}
			
			bool operator()() const { return enabled_; }

			size_t stackSize(size_t what) const
			{
				if (what==ENVIRON) return envStack_.size();
				return systemStack_.size();
			}

		private:
			const ParametersType& parameters_;
			bool enabled_;
			StackType systemStack_,envStack_; // <--we're the owner

			//! shrink  (we don't really shrink, we just undo the growth)
			void shrink(BasisType &pSprime,StackType& thisStack)
			{
				thisStack.pop();
				pSprime=thisStack.top();
			}
			
			//! Move elsewhere
			//! returns s1+s2 if s2 has no '/', 
			//! if s2 = s2a + '/' + s2b return s2a + '/' + s1 + s2b
			std::string appendWithDir(const std::string& s1,const std::string& s2) const
			{
				size_t x = s2.find("/");
				if (x==std::string::npos) return s1 + s2;
				std::string suf = s2.substr(x+1,s2.length());
				std::string dir = s2.substr(0,s2.length()-suf.length());
				//throw std::runtime_error("testing\n");
				return dir + s1 + suf;
			}	
	}; // class DmrgSerializer
} // namespace Dmrg 

/*@}*/
#endif
