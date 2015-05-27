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

/*! \file SymmetryElectronsSz.h
 *
 *  Stores info for each state of the Hilbert space basis
 *  This is a structure, don't add member functions!
 *
 */
#ifndef DMRG_SYMM_ELECTRONS_SZ_H
#define  DMRG_SYMM_ELECTRONS_SZ_H
#include "Vector.h"
#include "TypeToString.h"
#include "ProgramGlobals.h"
#include "IoSimple.h"
#include "TargetQuantumElectrons.h"

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
	         const VectorSizeType& szPlusConst)
	{
		jmValues_ = jmvalues;
		flavors_ = flavors;
		electrons_ = electrons;
		szPlusConst_ = szPlusConst;
	}

	SizeType electronsMax() const
	{
		return *(std::max_element(electrons_.begin(),electrons_.end()));
	}

	const VectorSizeType& electrons() const {return electrons_; }

	const VectorSizeType& flavors() const { return flavors_; }

	const VectorPairSizeType& jmValues() const
	{
		return jmValues_;
	}

	void findQuantumNumbers(VectorSizeType& qn, bool useSu2Symmetry) const
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
	                          const VectorSizeType& qns)
	{
		electrons.resize(qns.size());
		for (SizeType i=0;i<qns.size();i++) {
			VectorSizeType v = decodeQuantumNumber(qns[i]);
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

	static PsimagLite::String qnPrint(SizeType q)
	{
		PsimagLite::String str("");
		VectorSizeType qns = decodeQuantumNumber(q);
		for (SizeType k=0;k<qns.size();k++) str += ttos(qns[k]) + " ";
		return str;
	}

	static SizeType adjustQn(const VectorSizeType& adjustQuantumNumbers,
	                         SizeType direction,
	                         PsimagLite::IoSimple::Out& ioOut,
	                         bool useSu2Symmetry,
	                         SizeType step)
	{
		VectorSizeType targetQuantumNumbers(2,0);
		if (2*step+1 >= adjustQuantumNumbers.size()) {
			PsimagLite::String msg("adjustQuantumNumbers must be a vector");
			msg += " of size N-2, where N is the TotalNumberOfSites\n";
			throw PsimagLite::RuntimeError(msg);
		}

		targetQuantumNumbers[0] = adjustQuantumNumbers[2*step];
		targetQuantumNumbers[1] = adjustQuantumNumbers[2*step+1];
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

	static void setTargetNumbers(VectorSizeType& t,
	                             const TargetQuantumElectronsType& targetQ,
	                             SizeType sites,
	                             SizeType totalSites,
	                             SizeType direction)
	{
		t.resize((targetQ.isSu2) ? 3 : 2,0);

		if (direction == ProgramGlobals::INFINITE) {
			t[0] = static_cast<SizeType>(round(targetQ.szPlusConst*sites/totalSites));
			t[1] = static_cast<SizeType>(round(targetQ.totalElectrons*sites/totalSites));
		} else {
			t[0] = targetQ.szPlusConst;
			t[1] = targetQ.totalElectrons;
		}

		if (t.size() < 3) return;

		RealType jReal = targetQ.twiceJ*sites/static_cast<RealType>(totalSites);
		SizeType tmp = (direction == ProgramGlobals::INFINITE) ?
		            static_cast<SizeType>(round(jReal)) : targetQ.twiceJ;

		if (targetQ.totalElectrons%2==0) {
			if (tmp%2 != 0) tmp++;
		} else {
			if (tmp%2 == 0) tmp++;
		}

		t[2] = tmp;
	}

	void findQuantumNumbersSu2(VectorSizeType& q) const
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
	void findQuantumNumbersLocal(VectorSizeType& q) const
	{
		q.clear();
		VectorSizeType qn(2);
		for (SizeType i=0;i<electrons_.size();i++) {
			// n
			qn[1] = electrons_[i];
			// sz + const.
			qn[0] = szPlusConst_[i];

			q.push_back(encodeQuantumNumber(qn));
		}
	}

	static SizeType getQuantumSector(const VectorSizeType& targetQuantumNumbers,
	                                 SizeType direction,
	                                 PsimagLite::IoSimple::Out* ioOut,
	                                 bool useSu2Symmetry)
	{
		PsimagLite::OstringStream msg;
		msg<<"Integer target quantum numbers are: ";
		for (SizeType ii=0;ii<targetQuantumNumbers.size();ii++)
			msg<<targetQuantumNumbers[ii]<<" ";
		std::cout<<msg<<"\n";
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

		assert(v[0] < maxElectrons);
		assert(v[1] < maxElectrons);
		assert(v[2] < maxElectrons);

		SizeType x= v[0] + v[1]*maxElectrons;
		if (v.size()==3) x += v[2]*maxElectrons*maxElectrons;
		return x;
	}

	static VectorSizeType decodeQuantumNumber(SizeType q)
	{
		SizeType maxElectrons = 2*ProgramGlobals::maxElectronsOneSpin;

		assert(q < maxElectrons*maxElectrons*maxElectrons);

		VectorSizeType v(3);
		v[2] = SizeType(q/(maxElectrons*maxElectrons));
		SizeType tmp = q - v[2]*maxElectrons*maxElectrons;
		v[1] = SizeType(tmp/maxElectrons);
		v[0] = tmp % maxElectrons;
		return v;
	}

	VectorSizeType electrons_;
	VectorSizeType szPlusConst_;
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

