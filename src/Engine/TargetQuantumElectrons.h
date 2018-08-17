/*
Copyright (c) 2009-2015, UT-Battelle, LLC
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
template<typename RealType, typename QnType>
struct TargetQuantumElectrons {

	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename QnType::PairSizeType PairSizeType;

	template<typename IoInputType>
	TargetQuantumElectrons(IoInputType& io)
	    : totalNumberOfSites(0),
	      isSu2(false),
	      isCanonical(true),
	      qn(QnType::zero())
	{
		qn.other.clear();
		const bool allowUpDown = true;
		io.readline(totalNumberOfSites, "TotalNumberOfSites=");

		PsimagLite::String msg("TargetQuantumElectrons: ");
		bool hasTwiceJ = false;
		try {
			io.readline(qn.jmPair.first, "TargetSpinTimesTwo=");
			hasTwiceJ = true;
		} catch (std::exception&) {}

		SizeType ready = 0;
		if (allowUpDown) {
			SizeType electronsUp = 0;
			SizeType electronsDown = 0;
			try {
				io.readline(electronsUp,"TargetElectronsUp=");
				io.readline(electronsDown,"TargetElectronsDown=");
				SizeType tmp = electronsUp + electronsDown;
				qn.oddElectrons = (tmp & 1);
				qn.other.push_back(tmp);
				qn.other.push_back(electronsUp);
				ready = 2;
			} catch (std::exception&) {}
		}

		try {
			SizeType tmp = 0;
			io.readline(tmp, "TargetElectronsTotal=");
			qn.oddElectrons = (tmp & 1);
			qn.other.push_back(tmp);
			ready++;
		} catch (std::exception&) {}

		bool hasSzPlusConst = false;
		try {
			SizeType szPlusConst = 0;
			io.readline(szPlusConst,"TargetSzPlusConst=");
			qn.other.push_back(szPlusConst);
			hasSzPlusConst = true;
		} catch (std::exception&) {}

		if (ready == 3) {
			msg += "Provide either up/down or total/sz but not both.\n";
			throw PsimagLite::RuntimeError(msg);
		}

		if (qn.other.size() > 0) hasSzPlusConst = true;

		if (!hasSzPlusConst) {
			std::cerr<<"TargetQuantumElectrons: Grand Canonical\n";
			assert(qn.other.size() == 0);
			isCanonical = false;
		}

		bool flag = false;
		try {
			int dummy = 0;
			io.readline(dummy, "TargetExtra=");
			flag = true;
		} catch (std::exception&) {}

		if (flag) err("Instead of TargetExtra= please use a vector\n");

		try {
			VectorSizeType extra;
			io.read(extra,"TargetExtra");
			if (!hasSzPlusConst)
				std::cout<<"WARNING: TargetExtra= with grand canonical ???\n";
			for (SizeType i = 0; i < extra.size(); ++i)
				qn.other.push_back(extra[i]);
		} catch (std::exception&) {}

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

		if (isSu2)
			qn.oddElectrons = (totalNumberOfSites & 1);
	}

	template<typename SomeMemResolvType>
	SizeType memResolv(SomeMemResolvType&,
	                   SizeType,
	                   PsimagLite::String = "") const
	{
		return 0;
	}

	void write(PsimagLite::String label1,
	           PsimagLite::IoNg::Out::Serializer& io) const
	{
		PsimagLite::String label = label1 + "/TargetQuantumElectrons";
		io.createGroup(label);
		io.write(label + "/TotalNumberOfSites", totalNumberOfSites);
		io.write(label + "/isSu2", isSu2);
		io.write(label + "/isCanonical", isCanonical);
		qn.write(label + "/qn", io);
	}

	//! Function that prints model parameters to stream os
	friend std::ostream& operator<<(std::ostream &os,
	                                const TargetQuantumElectrons& p)
	{
		os<<"TargetElectronsTotal="<<p.totalElectrons<<"\n";
		os<<"TargetOther="<<p.other<<"\n";
		if (p.isSu2)
			os<<"TargetSpinTimesTwo="<<p.twiceJ<<"\n";
		return os;
	}

	SizeType totalNumberOfSites;
	bool isSu2;
	bool isCanonical;
	QnType qn;

private:

	TargetQuantumElectrons(const TargetQuantumElectrons&);

	TargetQuantumElectrons& operator=(const TargetQuantumElectrons&);
};
} // namespace Dmrg

/*@}*/
#endif

