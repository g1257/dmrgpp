/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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

/*! \file TargetQuantumElectrons.h
 *
 *
 */
#ifndef TargetQuantumElectrons_H
#define TargetQuantumElectrons_H
#include "Vector.h"
#include "ProgramGlobals.h"

namespace Dmrg {
//! Hubbard Model Parameters
template<typename RealType>
struct TargetQuantumElectrons {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	template<typename IoInputType>
	TargetQuantumElectrons(IoInputType& io, bool allowUpDown = true)
	    : twiceJ(0),isCanonical(true)
	{
		PsimagLite::String  msg("TargetQuantumElectrons: ");
		bool hasTwiceJ = false;
		try {
			io.readline(twiceJ,"TargetSpinTimesTwo=");
			hasTwiceJ = true;
		} catch (std::exception&) {}

		SizeType ready = 0;
		if (allowUpDown) {
			SizeType electronsUp = 0;
			SizeType electronsDown = 0;
			try {
				io.readline(electronsUp,"TargetElectronsUp=");
				io.readline(electronsDown,"TargetElectronsDown=");
				totalElectrons = electronsUp + electronsDown;
				other.push_back(electronsUp);
				ready=2;
			} catch (std::exception&) {}
		}

		try {
			io.readline(totalElectrons,"TargetElectronsTotal=");
			ready++;
		} catch (std::exception&) {}

		bool hasSzPlusConst = false;
		try {
			SizeType szPlusConst = 0;
			io.readline(szPlusConst,"TargetSzPlusConst=");
			other.push_back(szPlusConst);
			hasSzPlusConst = true;
		} catch (std::exception&) {}

		switch (ready) {
		case 3:
			msg += "Provide either up/down or total/sz but not both.\n";
			throw PsimagLite::RuntimeError(msg);

		case 0:
			msg += "Provide at least one of up/down or total/sz.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (other.size() > 0) hasSzPlusConst = true;

		if (!hasSzPlusConst) {
			std::cout<<"TargetQuantumElectrons: Grand Canonical\n";
			assert(other.size() == 0);
			other.resize(1,0);
			isCanonical = false;
		}

		while (true) {
			try {
				SizeType extra = 0;
				io.readline(extra,"TargetExtra=");
				if (other.size() == 0)
					std::cout<<"WARNING: TargetExtra= with grand canonical ???\n";
				other.push_back(extra);
			} catch (std::exception&) {
				break;
			}
		}

		int tmp = 0;
		try {
			io.readline(tmp,"UseSu2Symmetry=");
		} catch (std::exception&) {}

		isSu2 = (tmp > 0);

		if (isSu2 && !hasTwiceJ) {
			msg += "Please provide TargetSpinTimesTwo when running with SU(2).\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (isSu2 && !hasSzPlusConst)
			throw PsimagLite::RuntimeError
		        ("WARNING: SU(2) with grand canonical ???\n");
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	bool isSu2;
	SizeType totalElectrons;
	VectorSizeType other;
	SizeType twiceJ;
	bool isCanonical;
};

//! Function that prints model parameters to stream os
template<typename RealTypeType>
std::ostream& operator<<(std::ostream &os,
                         const TargetQuantumElectrons<RealTypeType>& p)
{
	os<<"TargetElectronsTotal="<<p.totalElectrons<<"\n";
	os<<"TargetOther="<<p.other<<"\n";
	if (p.isSu2)
		os<<"TargetSpinTimesTwo="<<p.twiceJ<<"\n";
	return os;
}
} // namespace Dmrg

/*@}*/
#endif

