/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file ParametersModelHubbard.h
 *
 *  Contains the parameters for the Hubbard model and function to read them from a file
 *
 */
#ifndef PARAMETERSMODELHUBBARD_H
#define PARAMETERSMODELHUBBARD_H

namespace Dmrg {
//! Hubbard Model Parameters
template<typename Field>
struct ParametersModelHubbard {

	template<typename IoInputType>
	ParametersModelHubbard(IoInputType& io)
	{

		io.read(hubbardU,"hubbardU");
		io.read(potentialV,"potentialV");
		//			SizeType level = 0;
		//			bool beQuiet = true;
		try {
			io.read(potentialT,"PotentialT"); //level,beQuiet);
		} catch (std::exception& e) {
		}
		omega=0;
		try {
			io.readline(omega,"omega=");
		} catch (std::exception& e) {
		}

		//			io.rewind();
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType& mres,
	                   SizeType,
	                   PsimagLite::String msg = "") const
	{
		PsimagLite::String str = msg;
		str += "ParametersModelHubbard";

		const char* start = reinterpret_cast<const char *>(this);
		const char* end = reinterpret_cast<const char *>(&potentialV);
		SizeType total = mres.memResolv(&hubbardU, end-start, str + " hubbardU");

		start = end;
		end = reinterpret_cast<const char *>(&potentialT);
		total += mres.memResolv(&potentialV, end-start, str + " potentialV");

		start = end;
		end = reinterpret_cast<const char *>(&omega);
		total += mres.memResolv(&potentialT, end-start, str + " potentialT");

		total += mres.memResolv(&omega,
		                        sizeof(*this) - total, str + " omega");

		return total;
	}

	//serializr start class ParametersModelHubbard
	// Do not include here connection parameters
	// those are handled by the Geometry
	// Hubbard U values (one for each site)
	//serializr normal hubbardU
	typename PsimagLite::Vector<Field>::Type hubbardU;
	// Onsite potential values, one for each site
	//serializr normal potentialV
	typename PsimagLite::Vector<Field>::Type potentialV;

	// for time-dependent H:
	//serializr normal potentialT
	typename PsimagLite::Vector<Field>::Type potentialT;
	//serializr normal omega
	Field omega;
};

//! Function that prints model parameters to stream os
template<typename FieldType>
std::ostream& operator<<(std::ostream &os,
                         const ParametersModelHubbard<FieldType>& parameters)
{
	os<<"hubbardU\n";
	os<<parameters.hubbardU;
	os<<"potentialV\n";
	os<<parameters.potentialV;
	if (parameters.potentialT.size()==0) return os;

	// time-dependent stuff
	os<<"potentialT\n";
	os<<parameters.potentialT;
	os<<"omega="<<parameters.omega<<"\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif

