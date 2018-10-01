/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

/*! \file ProgramGlobals.h
 *
 *
 *
 */
#ifndef PROGRAM_LIMITS_H
#define PROGRAM_LIMITS_H
#include "PsimagLite.h"
#include "Vector.h"
#include "Utils.h"
#include "../../src/Version.h"

namespace Dmrg {

struct ProgramGlobals {

	static SizeType maxElectronsOneSpin;

	static bool oldChangeOfBasis;

	static const PsimagLite::String license;

	static const SizeType MAX_LPS = 1000;

	static PsimagLite::String notReallySortAlgo;

	enum DirectionEnum {INFINITE=0,EXPAND_ENVIRON=1,EXPAND_SYSTEM=2};

	enum ConnectionEnum {SYSTEM_SYSTEM,SYSTEM_ENVIRON,ENVIRON_SYSTEM,ENVIRON_ENVIRON};

	enum SysOrEnvEnum {SYSTEM, ENVIRON};

	enum FermionOrBosonEnum {FERMION,BOSON};

	enum VerboseEnum {VERBOSE_NO, VERBOSE_YES};

	static void init(SizeType maxElectronsOneSpin_)
	{
		if (maxElectronsOneSpin == maxElectronsOneSpin_) return;
		if (maxElectronsOneSpin != 0)
			throw PsimagLite::RuntimeError("ProgramGlobals::init(...)\n");

		maxElectronsOneSpin = maxElectronsOneSpin_;
	}

	static int findBorderSiteFrom(SizeType site,
	                              DirectionEnum direction,
	                              SizeType n)
	{
		if (site == 1 && direction == EXPAND_ENVIRON)
			return 0;

		if (site == n - 2 && direction == EXPAND_SYSTEM)
			return n - 1;

		return -1;
	}

	static PsimagLite::String rootName(PsimagLite::String filename)
	{
		PsimagLite::String rootname = filename;
		size_t index =rootname.find(".", 0);
		if (index != PsimagLite::String::npos) {
			rootname.erase(index,filename.length());
		}

		return rootname;
	}

	static PsimagLite::String coutName(PsimagLite::String filename)
	{
		PsimagLite::String rootname = utils::basename(filename);
		size_t index =rootname.find(".", 0);
		if (index != PsimagLite::String::npos) {
			rootname.erase(index,filename.length());
		}

		return "runFor" + rootname + ".cout";
	}

	static SizeType logBase2(SizeType x)
	{
		SizeType counter = 0;
		while (x > 0) {
			x >>= 1;
			counter++;
		}

		return (counter == 0) ? counter : counter - 1;
	}

	static SizeType volumeOf(const PsimagLite::Vector<SizeType>::Type& v)
	{
		assert(v.size()>0);
		SizeType ret = v[0];
		for (SizeType i=1;i<v.size();i++) ret *= v[i];
		return ret;
	}

	friend std::istream& operator>>(std::istream& is, DirectionEnum& direction)
	{
		int x = -1;
		is >> x;
		if (x == 0) {
			direction = INFINITE;
		} else if (x == 1) {
			direction = EXPAND_ENVIRON;
		} else if (x == 2) {
			direction = EXPAND_SYSTEM;
		} else {
			err("istream& operator>> DirectionEnum\n");
		}

		return is;
	}

	static PsimagLite::String SYSTEM_STACK_STRING;
	static PsimagLite::String ENVIRON_STACK_STRING;
}; // ProgramGlobals

} // namespace Dmrg
/*@}*/
#endif

