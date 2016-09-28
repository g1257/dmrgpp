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

/*! \file SymmetryElectronsSz.h
 *
 *
 */
#ifndef DMRG_SYMM_ELECTRONS_SZ_H
#define  DMRG_SYMM_ELECTRONS_SZ_H
#include "Vector.h"
#include "TypeToString.h"
#include "ProgramGlobals.h"
#include "IoSimple.h"
#include "TargetQuantumElectrons.h"
#include "CvectorSize.h"

namespace Dmrg {

template<typename RealType>
class SymmetryElectronsSz {

	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairSizeType;
	typedef TargetQuantumElectrons<RealType> TargetQuantumElectronsType;

public:

	void set(const VectorPairSizeType& jmvalues,
	         const VectorSizeType& flavors,
	         const VectorSizeType& electrons,
	         const VectorSizeType& other)
	{
		jmValues_ = jmvalues;
		flavors_ = flavors;
		electrons_ = electrons;
		other_ = other;
	}

	SizeType electronsMax() const
	{
		return *(std::max_element(electrons_.begin(),electrons_.end()));
	}

	const CvectorSizeType& electrons() const {return electrons_; }

	const VectorSizeType& flavors() const { return flavors_; }

	const VectorPairSizeType& jmValues() const
	{
		return jmValues_;
	}

	void findQuantumNumbers(CvectorSizeType& qn, bool useSu2Symmetry) const
	{
		if (useSu2Symmetry)
			findQuantumNumbersSu2(qn);
		else
			findQuantumNumbersLocal(qn);
	}

	static SizeType neJmToIndex(SizeType ne,const PairType& jm)
	{
		VectorSizeType v(3);
		v[0] = jm.second;
		v[1] = ne;
		v[2] = jm.first;
		return encodeQuantumNumber(v);
	}

	static void qnToElectrons(VectorSizeType& electrons,
	                          const VectorSizeType& qns,
	                          SizeType total)
	{
		electrons.resize(qns.size());
		for (SizeType i=0;i<qns.size();i++) {
			VectorSizeType v = decodeQuantumNumber(qns[i],total);
			electrons[i] = v[1];
		}
	}

	static SizeType pseudoEffectiveNumber(SizeType nelectrons,
	                                      SizeType jtilde)
	{
		VectorSizeType v(3);
		v[0] = 0;
		v[1] = nelectrons;
		v[2] = jtilde;
		return pseudoQuantumNumber_(v);
	}

	static PsimagLite::String qnPrint(SizeType q, SizeType total)
	{
		PsimagLite::String str("");
		VectorSizeType qns = decodeQuantumNumber(q,total);
		for (SizeType k=0;k<qns.size();k++) str += ttos(qns[k]) + " ";
		return str;
	}

	static SizeType adjustQn(const VectorSizeType& adjustQuantumNumbers,
	                         SizeType direction,
	                         PsimagLite::IoSimple::Out& ioOut,
	                         bool useSu2Symmetry,
	                         SizeType step,
	                         SizeType mode)
	{
		VectorSizeType targetQuantumNumbers(mode+1,0);
		if (2*step+mode >= adjustQuantumNumbers.size()) {
			PsimagLite::String msg("adjustQuantumNumbers must be a vector");
			msg += " of correct size\n";
			throw PsimagLite::RuntimeError(msg);
		}

		for (SizeType x = 0; x < (mode+1); ++x)
			targetQuantumNumbers[x] = adjustQuantumNumbers[2*step+x];

		return getQuantumSector(targetQuantumNumbers,
		                        direction,
		                        &ioOut,
		                        useSu2Symmetry);
	}

	static SizeType getQuantumSector(const TargetQuantumElectronsType& targetQuantum,
	                                 SizeType sites,
	                                 SizeType total,
	                                 SizeType direction,
	                                 PsimagLite::IoSimple::Out* ioOut,
	                                 bool useSu2Symmetry)
	{
		VectorSizeType v;
		setTargetNumbers(v,targetQuantum,sites,total,direction);
		return getQuantumSector(v,direction,ioOut,useSu2Symmetry);
	}

private:

	void findQuantumNumbersSu2(CvectorSizeType& q) const
	{
		q.resize(electrons_.size());
		for (SizeType i=0;i<q.size();i++) {
			SizeType ne = electrons_[i];
			PairType jmpair = jmValues_[i];
			q[i]=neJmToIndex(ne,jmpair);
		}
	}

	//! find quantum numbers for each state of this basis,
	//! considered symmetries for this model are: n_up and n_down
	void findQuantumNumbersLocal(CvectorSizeType& q) const
	{
		SizeType mode = static_cast<SizeType>(other_.size()/electrons_.size());
		assert(other_.size() % electrons_.size() == 0);
		assert(mode > 0);

		q.clear();
		q.resize(electrons_.size());

		VectorSizeType qn(mode + 1);
		for (SizeType i=0;i<electrons_.size();i++) {
			// n
			qn[1] = electrons_[i];
			// sz + const.
			qn[0] = other_[i];

			for (SizeType x = 0; x < (mode-1); ++x)
				qn[2+x] = other_[i + (x+1)*electrons_.size()];
			q[i] = encodeQuantumNumber(qn);
		}
	}

	static void setTargetNumbers(VectorSizeType& t,
	                             const TargetQuantumElectronsType& targetQ,
	                             SizeType sites,
	                             SizeType totalSites,
	                             SizeType direction)
	{
		SizeType mode = targetQ.other.size();
		assert(mode > 0);
		assert(!targetQ.isSu2 || mode == 1);

		t.resize((targetQ.isSu2) ? 3 : mode + 1,0);

		if (direction == ProgramGlobals::INFINITE) {
			t[0] = static_cast<SizeType>(round(targetQ.other[0]*sites/totalSites));
			t[1] = static_cast<SizeType>(round(targetQ.totalElectrons*sites/totalSites));
			for (SizeType x = 0; x < (mode-1); ++x)
				t[2+x] = static_cast<SizeType>(round(targetQ.other[x+1]*sites/totalSites));
		} else {
			t[0] = targetQ.other[0];
			t[1] = targetQ.totalElectrons;
			for (SizeType x = 0; x < (mode-1); ++x)
				t[2+x] = static_cast<SizeType>(round(targetQ.other[x+1]*sites/totalSites));
		}

		if (!targetQ.isSu2) return;

		RealType jReal = targetQ.twiceJ*sites/static_cast<RealType>(totalSites);
		SizeType tmp = (direction == ProgramGlobals::INFINITE) ?
		            static_cast<SizeType>(round(jReal)) : targetQ.twiceJ;

		PsimagLite::String str("SymmetryElectronsSz: FATAL: Impossible parameters ");
		bool flag = false;
		if (targetQ.totalElectrons & 1) {
			if (!(tmp&1)) {
				flag = true;
				str += "electrons= " + ttos(targetQ.totalElectrons) + " is odd ";
				str += "and 2j= " +  ttos(tmp) + " is even.";
				tmp++;
			}
		} else {
			if (tmp & 1) {
				flag = true;
				str += "electrons= " + ttos(targetQ.totalElectrons) + " is even ";
				str += "and 2j= " +  ttos(tmp) + " is odd.";
				tmp++;
			}
		}

		if (flag && sites == totalSites) throw PsimagLite::RuntimeError(str);

		t[2] = tmp;
	}

	static SizeType getQuantumSector(const VectorSizeType& targetQuantumNumbers,
	                                 SizeType direction,
	                                 PsimagLite::IoSimple::Out* ioOut,
	                                 bool useSu2Symmetry)
	{
		PsimagLite::OstringStream msg;
		msg<<"SymmetryElectronsSz: Integer target quantum numbers are: ";
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		std::cout<<msg.str()<<"\n";
		if (ioOut && direction == ProgramGlobals::INFINITE)
			ioOut->printVector(targetQuantumNumbers,"TargetedQuantumNumbers");
		return pseudoQuantumNumber(targetQuantumNumbers,useSu2Symmetry);
	}

	//! Encodes (flavor,jvalue,density) into a unique number and returns it
	static SizeType pseudoQuantumNumber(const VectorSizeType& targets,
	                                    bool useSu2Symmetry)
	{
		if (useSu2Symmetry)
			return pseudoQuantumNumber_(targets);
		else
			return encodeQuantumNumber(targets);
	}

	//! targets[0]=nup, targets[1]=ndown,  targets[2]=2j
	static SizeType pseudoQuantumNumber_(const VectorSizeType& v)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		SizeType x = v[1];

		assert(x < maxElectrons);

		x += v[2]*maxElectrons;
		return x;
	}

	static SizeType encodeQuantumNumber(const VectorSizeType& v)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		for (SizeType x = 0; x < v.size(); ++x)
			assert(v[x] < maxElectrons);

		SizeType index = 0;
		SizeType number = 1;
		for (SizeType x = 0; x < v.size(); ++x) {
			index += v[x] * number;
			number *= maxElectrons;
		}

		return index;
	}

	static VectorSizeType decodeQuantumNumber(SizeType q, SizeType total)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		SizeType number = 1;
		for (SizeType x = 1; x < total; ++x) number *= maxElectrons;
		SizeType tmp = q;
		VectorSizeType v(total);
		for (SizeType x = 0; x < total; ++x) {
			v[total-1-x] = static_cast<SizeType>(tmp/number);
			tmp -= v[total-1-x]*number;
			number /= maxElectrons;

		}

		return v;
	}

	CvectorSizeType electrons_;
	VectorSizeType other_;
	VectorPairSizeType jmValues_;
	VectorSizeType flavors_;
}; // struct SymmetryElectronsSz

template<typename IoOutputter, typename RealType>
void save(IoOutputter& io,const SymmetryElectronsSz<RealType>& bd)
{
	io.printVector(bd.electronsUp,"#bdElectronsUp");
	io.printVector(bd.electronsDown,"#bdElectronsDown");
	io.printVector(bd.jmValues,"#bdJmValues");
	io.printVector(bd.flavors,"#bdFlavors=");
	PsimagLite::String msg("Symmetry=ElectronsSz\n");
	io.print(msg);
}

} // namespace Dmrg

/*@}*/
#endif

