/*
Copyright (c) 2009-2016, UT-Battelle, LLC
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

/*! \file Su2Related.h
 *
 *
 */
#ifndef DMRG_SU2_RELATED_H
#define DMRG_SU2_RELATED_H
#include "CrsMatrix.h"
#include "InputNg.h"
#include "InputCheck.h"

namespace Dmrg {
// This is a structure, don't add member functions here!
struct Su2Related {
	Su2Related()
	    : offset(0) // setting to zero is necessary, because
	{}		// we always print offset
	// and when running Abelian
	// it might be undefined

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "Su2Related";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&source);
		SizeType total = mres.memResolv(&offset, end-start, str + " offset");

		start = end;
		end = reinterpret_cast<const char *>(&transpose);
		total += mres.memResolv(&source, end-start, str + " source");

		total += mres.memResolv(&transpose,
		                        sizeof(*this) - total,
		                        str + " transpose");

		return total;
	}

	//serializr start class Su2Related
	//serializr normal offset
	SizeType offset;
	//serializr normal source
	PsimagLite::Vector<SizeType>::Type source;
	//serializr normal transpose
	PsimagLite::Vector<int>::Type transpose;
};

std::istream& operator>>(std::istream&,Su2Related&);

std::ostream& operator<<(std::ostream&,const Su2Related&);

void send(Su2Related&,int,int,PsimagLite::MPI::CommType);

void recv(Su2Related&,int,int,PsimagLite::MPI::CommType);

void bcast(Su2Related&);

} // namespace Dmrg

#endif

