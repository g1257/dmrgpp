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

/*! \file LinkProductTjMultiOrb.h
 *
 *  FIXME
 *
 */
#ifndef DMRG_LINK_PROD_TJ_MULTIORB_H
#define DMRG_LINK_PROD_TJ_MULTIORB_H


namespace Dmrg {

template<typename ModelHelperType>
class LinkProductTjMultiOrb {
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef std::pair<SizeType,SizeType> PairType;

	enum {TERM_CICJ, TERM_SISJ, TERM_NINJ};

public:
	typedef typename ModelHelperType::RealType RealType;
	typedef typename SparseMatrixType::value_type SparseElementType;

	template<typename SomeStructType>
	static void setLinkData(SizeType term,
							SizeType dofs,
							bool isSu2,
							SizeType& fermionOrBoson,
							PairType& ops,
							std::pair<char,char>& mods,
							SizeType& angularMomentum,
							RealType& angularFactor,
							SizeType& category,
							const SomeStructType&)
	{
		char tmp = mods.first;
		if (term==TERM_CICJ) {
			fermionOrBoson = ProgramGlobals::FERMION;
			switch (dofs) {
			case 0: // up up
				angularFactor = 1;
				category = 0;
				angularMomentum = 1;
				ops = PairType(dofs,dofs);
				break;
			case 1: // down down
				angularFactor = -1;
				category = 1;
				angularMomentum = 1;
				ops = PairType(dofs,dofs);
				break;
			}
			return;
		}

		if (term==TERM_SISJ) {
			fermionOrBoson = ProgramGlobals::BOSON;
			switch (dofs) {
			case 0: // S+ S-
				angularFactor = -1;
				category = 2;
				angularMomentum = 2;
				ops = operatorDofs(dofs,isSu2);
				break;
			case 1: // S- S+
				angularFactor = -1;
				category = 0;
				mods.first = mods.second;
				mods.second = tmp;
				angularMomentum = 2;
				ops = operatorDofs(dofs,isSu2);
				break;
			case 2: // Sz Sz
				angularFactor = 0.5;
				category = 1;
				angularMomentum = 2;
				ops = operatorDofs(dofs,isSu2);
				break;
			}
			return;
		}

		if (term==TERM_NINJ) {
			fermionOrBoson = ProgramGlobals::BOSON;
			ops = PairType(4,4);
			angularFactor = 1;
			angularMomentum = 0;
			category = 0;
			return;
		}

		assert(false);
	}

	template<typename SomeStructType>
	static void valueModifier(SparseElementType& value,SizeType term,SizeType dofs,bool isSu2,const SomeStructType&)
	{
		if (term==TERM_CICJ) return;

		if (term==TERM_NINJ) {
			value *= 0.5;
			return;
		}

		assert(term==TERM_SISJ);

		if (isSu2) value = -value;
		if (dofs<2) value *= 0.5;
		value *= 0.5;
	}

	// connections are:
	// up up and down down
	// S+ S- and S- S+
	// Sz Sz
	template<typename SomeStructType>
	static SizeType dofs(SizeType term,const SomeStructType&)
	{
		if (term==TERM_CICJ) return 2; // c^\dagger c
		if (term==TERM_SISJ) return 3; // S+ S-, S- S+ and Sz Sz
		if (term==TERM_NINJ) return 1; // ninj
		assert(false);
		return 0; // bogus
	}

	template<typename SomeStructType>
	static std::pair<SizeType,SizeType> connectorDofs(SizeType,SizeType,const SomeStructType&)
	{
		return PairType(0,0); // no orbital and no dependence on spin
	}

	static SizeType terms() { return 2; }

private:

	// only for TERM_SISJ
	static PairType operatorDofs(SizeType dofs,bool isSu2)
	{
		if (dofs<2) return PairType(2,2);
		SizeType x = (isSu2) ? 2 : 3;
		return PairType(x,x);
	}

}; // class LinkProductTjMultiOrb
} // namespace Dmrg
/*@}*/
#endif
