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

/*! \file LinkProductHeisenbergAncillaC.h
 *
 *  LinkProduct for Heisenberg model
 *
 */
#ifndef DMRG_LINK_PROD_HEISENBERG_ANCILLAC_H
#define DMRG_LINK_PROD_HEISENBERG_ANCILLAC_H
#include "ProgramGlobals.h"

namespace Dmrg {
template<typename ModelHelperType>
class LinkProductHeisenbergAncillaC {

public:

	enum {TERM_SPLUSSMINUS, TERM_SZSZ, TERM_ANCILLA};

	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename ModelHelperType::RealType RealType;

	template<typename SomeStructType>
	static void setLinkData(SizeType term,
	                        SizeType dofs,
	                        bool isSu2,
	                        ProgramGlobals::FermionOrBosonEnum& fermionOrBoson,
	                        std::pair<SizeType,SizeType>& ops,
	                        std::pair<char,char>&,
	                        SizeType& angularMomentum,
	                        RealType& angularFactor,
	                        SizeType& category,
	                        const SomeStructType&)
	{
		fermionOrBoson = ProgramGlobals::BOSON;
		ops = (hot_) ? operatorDofsHot(term,dofs,isSu2) : operatorDofs(term,isSu2);
		angularMomentum = 2;
		switch (term) {
		case 0:
			angularFactor = -1;
			category = 2;
			break;
		case 1:
			angularFactor = 0.5;
			category = 1;
			break;
		case 2:
			angularFactor = 1;
			category = 0;
			break;
		}
	}

	template<typename SomeStructType>
	static void valueModifier(SparseElementType& value,
	                          SizeType term,
	                          SizeType,
	                          bool isSu2,
	                          const SomeStructType&)
	{
		if (term == TERM_ANCILLA) return;
		if (isSu2) value = -value;
		value *= 0.5;
	}

	template<typename SomeStructType>
	static SizeType dofs(SizeType term,const SomeStructType&)
	{
		if (!hot_ || term == TERM_ANCILLA) return 1;
		return 2;
	}

	template<typename SomeStructType>
	static PairType connectorDofs(SizeType term,SizeType dofs,const SomeStructType&)
	{
		if (!hot_ || term == TERM_ANCILLA) return PairType(0,0);
		return PairType(dofs,dofs);
	}

	//! Splus Sminus and
	//! Sz Sz
	//! delta^\dagger delta
	static SizeType terms() { return 3; }

	static bool setHot(bool hot) { return hot_ = hot; }

private:

	static PairType operatorDofs(SizeType term,bool isSu2)
	{
		if (term == TERM_SPLUSSMINUS || term == TERM_ANCILLA)
			return PairType(term,term);
		SizeType x = (isSu2) ? 0 : 1;
		return PairType(x,x);
	}

	static PairType operatorDofsHot(SizeType term,SizeType dofs,bool)
	{
		if (term == TERM_ANCILLA) return PairType(4,4);
		assert(term == TERM_SPLUSSMINUS || term == TERM_SZSZ);
		SizeType offset = (term == TERM_SPLUSSMINUS) ? 0 : 2;
		return PairType(dofs+offset,dofs+offset);
	}

	static bool hot_;
}; // class LinkProductHeisenbergAncillaC
} // namespace Dmrg
/*@}*/
#endif

