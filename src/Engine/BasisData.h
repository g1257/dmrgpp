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

/*! \file BasisData.h
 *
 *  Stores info for each state of the Hilbert space basis
 *  This is a structure, don't add member functions!
 *
 */
#ifndef BASIS_DATA_H
#define  BASIS_DATA_H
#include "Vector.h"


namespace Dmrg {
template<typename PairType>
class BasisData {

	typedef typename PsimagLite::Vector<SizeType>::Type VectorSizeType;
	typedef typename PsimagLite::Vector<PairType>::Type VectorPairSizeType;

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

			//assert(qn[1]>=qn[0]);
			//qn[1] -= qn[0];

			q.push_back(encodeQuantumNumber(qn));
		}
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

	static SizeType neJmToIndex(SizeType ne,const PairType& jm)
	{
		VectorSizeType v(3);
		v[0] = jm.second;
		v[1] = ne;
		v[2] = jm.first;
		return encodeQuantumNumber(v);
	}

private:

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

	VectorSizeType electrons_;
	VectorSizeType szPlusConst_;
	VectorPairSizeType jmValues_;
	VectorSizeType flavors_;
}; // struct BasisData

template<typename IoOutputter, typename PairType>
void save(IoOutputter& io,const BasisData<PairType>& bd)
{
	io.printVector(bd.electronsUp,"#bdElectronsUp");
	io.printVector(bd.electronsDown,"#bdElectronsDown");
	io.printVector(bd.jmValues,"#bdJmValues");
	io.printVector(bd.flavors,"#bdFlavors=");
	PsimagLite::String msg("Symmetry=");
	msg += bd.symm->name();
}

} // namespace Dmrg

/*@}*/
#endif

