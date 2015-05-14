/*
Copyright (c) 2009-2014, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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

	bool check(const PsimagLite::String& label,
	           const PsimagLite::Vector<PsimagLite::String>::Type& vec,
	           SizeType line) const
	{
		if (label=="JMVALUES") {
			if (vec.size()!=2) return error1("JMVALUES",line);
			return true;
		} else if (label=="RAW_MATRIX") {
			if (!checkForMatrix(vec)) return error1(label,line);
			return true;
		} else if (label=="Connectors") {
			if (!checkForMatrix(vec) && !checkForVector(vec))
				return error1(label,line);
			return true;
		} else if (label=="MagneticField") {
			return true;
		} else if (label=="FiniteLoops") {
			SizeType n = atoi(vec[0].c_str());
			if (vec.size()!=3*n+1)  return error1("FiniteLoops",line);
			return true;
		}
		return false;
	}

	/* PSIDOC dmrgSolverOptions
		  \item[Options=string]
		  A comma-separated list of strings. At least one of the following strings must
		  be provided:
		  \begin{itemize}
			 \item[none]  Use this when no options are given, since the list of
		   strings must be non-null.
				Note that ``none'' does not disable other options.

			 \item[useSu2Symmetry] Use the SU(2) symmetry for the model, and
			interpret quantum
				 numbers in the line ``QNS'' appropriately.

			 \item[nofiniteloops]  Don't do finite loops, even if provided under
			``FiniteLoops'' below.
			\end{itemize}
		*/
	void check(const PsimagLite::String& label,const PsimagLite::String& val,SizeType)
	{
		if (label!="SolverOptions") return;
		PsimagLite::Vector<PsimagLite::String>::Type registerOpts;

		registerOpts.push_back("restart");
		registerOpts.push_back("debugmatrix");
		registerOpts.push_back("test");
		registerOpts.push_back("exactdiag");
		registerOpts.push_back("useDavidson");
		registerOpts.push_back("verbose");
		registerOpts.push_back("nofiniteloops");
		registerOpts.push_back("nowft");
		registerOpts.push_back("targetnoguess");
		registerOpts.push_back("complex");
		registerOpts.push_back("inflate");
		registerOpts.push_back("none");
		registerOpts.push_back("twositedmrg");
		registerOpts.push_back("noloadwft");
		registerOpts.push_back("concurrenttridiag");
		registerOpts.push_back("ChebyshevSolver");
		registerOpts.push_back("MatrixVectorStored");
		registerOpts.push_back("MatrixVectorKron");
		registerOpts.push_back("TimeStepTargetting");
		registerOpts.push_back("DynamicTargetting");
		registerOpts.push_back("AdaptiveDynamicTargetting");
		registerOpts.push_back("CorrectionVectorTargetting");
		registerOpts.push_back("CorrectionTargetting");
		registerOpts.push_back("MettsTargetting");
		registerOpts.push_back("TargetingAncilla");
		registerOpts.push_back("geometryallinsystem");

		PsimagLite::Options::Writeable
		        optWriteable(registerOpts,PsimagLite::Options::Writeable::PERMISSIVE);
		optsReadable_ = new  OptionsReadableType(optWriteable,val);
	}

	bool isSet(const PsimagLite::String& thisOption) const
	{
		return optsReadable_->isSet(thisOption);
	}

	void checkForThreads(SizeType nthreads) const
	{
		if (nthreads==1) return;

		PsimagLite::String message1(__FILE__);
		message1 += " FATAL: You are requesting nthreads>0 but you ";
		message1 += "did not compile with USE_PTHREADS enabled\n";
		message1 += " Either set Threads=1 in the input file (you won't ";
		message1 += "have threads though) or\n";
		message1 += " add -DUSE_PTHREADS to the CPP_FLAGS in your Makefile ";
		message1 += "and recompile\n";
		throw PsimagLite::RuntimeError(message1.c_str());
	}

	void usageMain(const PsimagLite::String& name) const
	{
		std::cerr<<"USAGE is "<<name<<"\n";
	}

	PsimagLite::String getTargeting(const PsimagLite::String& options) const
	{
		PsimagLite::String targetting="GroundStateTargetting";

		const char *targets[]={"GroundStateTargetting",
		                       "TimeStepTargetting",
		                       "AdaptiveDynamicTargetting",
		                       "DynamicTargetting",
		                       "CorrectionVectorTargetting",
		                       "CorrectionTargetting",
		                       "MettsTargetting",
		                       "TargetingAncilla"};

		SizeType totalTargets = 8;

		SizeType count = 0;
		for (SizeType i = 0;i<totalTargets;++i) {
			if (options.find(targets[i])!=PsimagLite::String::npos) {
				if (targetting == "AdaptiveDynamicTargetting" &&
				        std::string(targets[i]) == "DynamicTargetting") continue;
				targetting = targets[i];
				count++;
			}
		}

		if (count > 1) {
			throw PsimagLite::RuntimeError("Only one targeting supported\n");
		}

		return targetting;
	}

private:

	bool checkForVector(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() == 0) return false;
		SizeType n = atoi(vec[0].c_str());
		return (vec.size() == n+1);
	}

	bool checkForMatrix(const PsimagLite::Vector<PsimagLite::String>::Type& vec) const
	{
		if (vec.size() < 2) return false;
		SizeType row = atoi(vec[0].c_str());
		SizeType col = atoi(vec[1].c_str());
		SizeType n = row*col;
		return (vec.size() == n+2);
	}

	bool error1(const PsimagLite::String& message,SizeType line) const
	{
		PsimagLite::String s(__FILE__);
		s += " : Input error for label " + message + " near line " + ttos(line) + "\n";
		throw PsimagLite::RuntimeError(s.c_str());

	}

	OptionsReadableType* optsReadable_;

}; // class InputCheck
} // namespace Dmrg

/*@}*/
#endif

