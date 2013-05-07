/*
Copyright (c) 2012 , UT-Battelle, LLC
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

/*! \file Options.h
 *
 *  Options lines
 */

#ifndef OPTIONS_HEADER_H
#define OPTIONS_HEADER_H

#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <iostream>

namespace PsimagLite {

class Options {

public:
	class Writeable {

	public:

		enum {DISABLED,PERMISSIVE,STRICT};
	
		Writeable(Vector<std::string>::Type& registeredOptions,size_t mode)
		: registeredOptions_(registeredOptions),mode_(mode)
		{}

		void set(Vector<std::string>::Type& optsThatAreSet,const std::string& opts)
		{
			if (mode_==DISABLED) return;
			split(optsThatAreSet,opts.c_str(),',');
			for (size_t i=0;i<optsThatAreSet.size();i++) {
				bool b = (find(registeredOptions_.begin(),registeredOptions_.end(),optsThatAreSet[i])==registeredOptions_.end());
				if (!b) continue;
				
				std::string s(__FILE__);
				s += ": Unknown option " + optsThatAreSet[i] + "\n";
				if (mode_==PERMISSIVE) std::cout<<" *** WARNING **: "<<s;
				if (mode_==STRICT) throw std::runtime_error(s.c_str());
			}
		}

	private:
		Vector<std::string>::Type registeredOptions_;
		size_t mode_;
	}; // class Writeable

	class Readable {

	public:
		Readable(Writeable& optsWrite,const std::string& optsString)
		{
			optsWrite.set(optsThatAreSet_,optsString);
		}

		bool isSet(const std::string& thisOption) const
		{
			bool b = (find(optsThatAreSet_.begin(),optsThatAreSet_.end(),thisOption)==optsThatAreSet_.end());
			return (!b);
		}
	private:
		Vector<std::string>::Type optsThatAreSet_;

	}; // class Readable

}; // class OptionsWriteable
} // namespace PsimagLite
/*@}*/
#endif // OPTIONS_HEADER_H
