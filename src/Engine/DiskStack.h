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
#ifndef DISKSTACK_HEADER_H
#define DISKSTACK_HEADER_H

// All these includes are in PsimagLite
#include "Stack.h"
#include "IoSimple.h"
#include "ProgressIndicator.h"

//! A disk stack, similar to std::stack but stores in disk not in memory
namespace Dmrg {
	template<typename DataType>
	class DiskStack {
	
		typedef typename PsimagLite::IoSimple::In IoInType;
		typedef typename PsimagLite::IoSimple::Out IoOutType;

		public:
			DiskStack(const std::string &file1,const std::string &file2,bool hasLoad,size_t rank=0) :
				rank_(rank),
				fileIn_(file1),
				fileOut_(file2),
				total_(0),
				progress_("DiskStack",rank)
			{
				if (!hasLoad) {
					ioOut_.open(fileOut_,std::ios_base::trunc,rank_);
					ioOut_.close();
					return;
				}
				try {
					ioIn_.open(fileIn_);
				} catch (std::exception& e) {
					std::cerr<<"Problem opening reading file "<<fileIn_<<"\n";
					throw std::runtime_error("DiskStack::load(...)\n");
				}
				ioIn_.readline(rank_,"#STACKMETARANK=",IoInType::LAST_INSTANCE);
				//ioIn_.readline(total_,"#STACKMETATOTAL="); <-- KEEP TOTAL=0, THIS IS A NEW FILE!!
				//ioIn_.readline(debug_,"#STACKMETADEBUG=");
				ioIn_.advance("#STACKMETASTACK");
				ioIn_>>stack_;
				ioIn_.close();
				std::ostringstream msg;
				msg<<"Attempting to read from file " + fileIn_ + " succeeded";
				progress_.printline(msg,std::cout);
				
			}
			
			~DiskStack()
			{
				//ioOut_.open(fileOut_,std::ios_base::trunc,rank_);
				ioOut_.open(fileOut_,std::ios_base::app,rank_);
				ioOut_.print("#STACKMETARANK=",rank_);
				
				ioOut_.print("#STACKMETATOTAL=",total_);
				//ioOut_.printline("#STACKMETADEBUG="+ttos(debug_));
				ioOut_<<"#STACKMETASTACK\n";
				ioOut_<<stack_;
				ioOut_.close();
			}

			static bool persistent() { return true; }

			void push(DataType const &d) 
			{
				//std::string tmpLabel = fileOut_ + ttos(total_);
				ioOut_.open(fileOut_,std::ios_base::app,rank_);
				d.save(ioOut_);
				ioOut_.close();

				stack_.push(total_);
				total_++;

				//std::string s = "Pushing with label="+fileOut_+" and total="+ttos(total_);
				//debugPrint(s);
				//std::ostringstream msg;
				//msg<<s;
				//progress_.printline(msg,std::cout);
			}

			void pop()
			{
				stack_.pop();
//				std::string s = "Popping with label="+fileIn_+" stack_.top="+ttos(stack_.top());
//				std::ostringstream msg;
//				msg<<s;
//				progress_.printline(msg,std::cout);
			}

			DataType top()
			{
				ioIn_.open(fileIn_);
				DataType dt(ioIn_,"",stack_.top());
//				std::string s = "Topping with label="+fileIn_+" stack_.top="+ttos(stack_.top());
//				std::ostringstream msg;
//				msg<<s;
//				progress_.printline(msg,std::cout);
				ioIn_.close();
				return dt;
			}

			size_t size() const { return stack_.size(); }

			template<typename DataType_>
			friend std::ostream& operator<<(std::ostream& os,const DiskStack<DataType_>& ds);

		private:

			size_t rank_;
			std::string fileIn_,fileOut_;
			int total_;
			PsimagLite::ProgressIndicator progress_;
			IoInType ioIn_;
			IoOutType ioOut_;
			std::stack<int> stack_;
	}; // class DiskStack

	template<typename DataType>
	std::ostream& operator<<(std::ostream& os,const DiskStack<DataType>& ds)
	{
		os<<"DISKSTACK: filein: "<<ds.fileIn_<<" fileout="<<ds.fileOut_<<"\n";
		os<<"total="<<ds.total_<<"\n";
		os<<ds.stack_;
		return os;
	}
} // namespace DMrg

#endif
