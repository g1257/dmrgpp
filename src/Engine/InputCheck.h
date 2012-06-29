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
/** \ingroup DMRG */
/*@{*/

/*! \file InputCheck.h
 *
 *  InputChecking functions
 */
#ifndef INPUT_CHECK_H
#define INPUT_CHECK_H
#include <vector>
#include <string>
#include <stdexcept>
#include "Options.h"

namespace Dmrg {

	class InputCheck {

		typedef PsimagLite::Options::Readable OptionsReadableType;

	public:

		InputCheck() : optsReadable_(0) {}

		~InputCheck()
		{
			if (optsReadable_!=0) delete optsReadable_;
		}

		bool check(const std::string& label,const std::vector<std::string>& vec,size_t line) const
		{
			if (label=="JMVALUES") {
				if (vec.size()!=2) return error1("JMVALUES",line);
				return true;
			} else if (label=="RAW_MATRIX") {
				size_t row = atoi(vec[0].c_str());
				size_t col = atoi(vec[1].c_str());
				size_t n = row*col;
				if (vec.size()!=n+2) return error1("RAW_MATRIX",line);
				return true;
			} else if (label=="Connectors") {
				return true;
			} else if (label=="MagneticField") {
				return true;
			} else if (label=="FiniteLoops") {
				size_t n = atoi(vec[0].c_str());
				if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
				return true;
			}
			return false;
		}

		void check(const std::string& label,const std::string& val,size_t line)
		{
			if (label!="SolverOptions") return;
			std::vector<std::string> registerOpts;

			registerOpts.push_back("restart");
			registerOpts.push_back("debugmatrix");
			registerOpts.push_back("test");
			registerOpts.push_back("useDavidson");
			registerOpts.push_back("verbose");
			registerOpts.push_back("nofiniteloops");
			registerOpts.push_back("nowft");
			registerOpts.push_back("inflate");
			registerOpts.push_back("none");
			registerOpts.push_back("ChebyshevSolver");
			registerOpts.push_back("InternalProductStored");
			registerOpts.push_back("InternalProductKron");
			registerOpts.push_back("useSu2Symmetry");
			registerOpts.push_back("TimeStepTargetting");
			registerOpts.push_back("DynamicTargetting");
			registerOpts.push_back("AdaptiveDynamicTargetting");
			registerOpts.push_back("CorrectionVectorTargetting");
			registerOpts.push_back("CorrectionTargetting");
			registerOpts.push_back("MettsTargetting");

			PsimagLite::Options::Writeable optWriteable(registerOpts,PsimagLite::Options::Writeable::PERMISSIVE);
			optsReadable_ = new  OptionsReadableType(optWriteable,val);
		}

		bool isSet(const std::string& thisOption) const
		{
			return optsReadable_->isSet(thisOption);
		}

		void checkForThreads(size_t nthreads) const
		{
			if (nthreads==1) return;

			std::string message1(__FILE__);
			message1 += " FATAL: You are requesting nthreads>0 but you did not compile with USE_PTHREADS enabled\n";
			message1 += " Either set Threads=1 in the input file (you won't have threads though) or\n";
			message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile and recompile\n";
			throw std::runtime_error(message1.c_str());
		}

		void usageMain(const std::string& name) const
		{
			std::cerr<<"USAGE is "<<name<<"\n";
		}

	private:

		bool error1(const std::string& message,size_t line) const
		{
			std::string s(__FILE__);
			s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
			throw std::runtime_error(s.c_str());

		}

		OptionsReadableType* optsReadable_;

	}; // class InputCheck
} // namespace Dmrg 

/*@}*/
#endif
