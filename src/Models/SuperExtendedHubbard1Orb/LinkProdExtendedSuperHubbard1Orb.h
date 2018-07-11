/*
Copyright (c) 2009-2016-2018, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file LinkProdExtendedSuperHubbard1Orb.h
 *
 *  FIXME
 *
 */
#ifndef LINKPROD_EXTENDED_SUPER_HUBBARD_1ORB_H
#define LINKPROD_EXTENDED_SUPER_HUBBARD_1ORB_H
#include "ProgramGlobals.h"
#include "LinkProductBase.h"

namespace Dmrg {

template<typename ModelHelperType, typename GeometryType>
class LinkProdExtendedSuperHubbard1Orb : public LinkProductBase<ModelHelperType, GeometryType> {

	typedef LinkProductBase<ModelHelperType, GeometryType> BaseType;
	typedef typename BaseType::AdditionalDataType AdditionalDataType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef std::pair<SizeType,SizeType> PairType;

public:

	enum {TERM_HOPPING=0,TERM_NINJ=1,TERM_SPLUSISMINUSJ=2,TERM_SZISZJ=3,TERM_PAIRIPAIRJ=4};

	typedef typename ModelHelperType::RealType RealType;
	typedef typename SparseMatrixType::value_type SparseElementType;

	void setLinkData(SizeType term,
	                 SizeType dofs,
	                 bool,
	                 ProgramGlobals::FermionOrBosonEnum& fermionOrBoson,
	                 PairType& ops,
	                 std::pair<char,char>&,
	                 SizeType& angularMomentum,
	                 RealType& angularFactor,
	                 SizeType& category,
	                 const AdditionalDataType&) const
	{
		if (term==TERM_HOPPING) {
			fermionOrBoson = ProgramGlobals::FERMION;
			SizeType spin1 = (dofs & 1);
			SizeType spin2 = (dofs & 2);
			spin2 >>= 1;
			if (!BaseType::hasSpinOrbit())
				spin1 = spin2 = dofs;
			ops = PairType(spin1,spin2);
			angularFactor = 1;
			if (spin1 == 1) angularFactor = -1;
			angularMomentum = 1;
			category = spin1;
			return;
		}

		if (term==TERM_NINJ) {
			fermionOrBoson = ProgramGlobals::BOSON;
			SizeType offset1 = 2;
			assert(dofs == 0);
			angularFactor = 1;
			category = 0;
			angularMomentum = 0;
			ops = PairType(offset1,offset1);
			return;
		}

		if (term==TERM_SPLUSISMINUSJ) {
			fermionOrBoson = ProgramGlobals::BOSON;
			SizeType offset1 = 3;
			angularFactor = -1;
			category = 2;
			angularMomentum = 2;
			ops = PairType(offset1,offset1);
			return;
		}

		if (term==TERM_SZISZJ) {
			fermionOrBoson = ProgramGlobals::BOSON;
			SizeType offset1 = 4;
			assert(dofs == 0);
			angularFactor = 0.5;
			category = 1;
			angularMomentum = 2;
			ops = PairType(offset1,offset1);
			return;
		}

		if (term==TERM_PAIRIPAIRJ) {
			fermionOrBoson = ProgramGlobals::BOSON;
			SizeType offset1 = 5;
			angularFactor = 1;
			angularMomentum = 2;
			category = 2;
			ops = PairType(offset1,offset1);
			return;
		}
	}

	void valueModifier(SparseElementType& value,
	                   SizeType term,
	                   SizeType,
	                   bool isSu2,
	                   const AdditionalDataType&) const
	{

		if (term==TERM_HOPPING || term == TERM_PAIRIPAIRJ)
			return;

		if (term==TERM_NINJ) {
			value *= 0.5;
			return;
		}

		assert(term==TERM_SPLUSISMINUSJ || term == TERM_SZISZJ);

		if (isSu2) value = -value;
		value *= 0.5;
	}

	SizeType dofs(SizeType term,const AdditionalDataType&) const
	{
		if (term != TERM_HOPPING) return 1;
		assert(term == TERM_HOPPING);
		return (BaseType::hasSpinOrbit()) ? 4 : 2;
	}

	// has only dependence on orbital
	void connectorDofs(VectorSizeType& edofs,
	                   SizeType term,
	                   SizeType dofs,
	                   const AdditionalDataType&) const
	{
		// no orbital and no dependence on spin
		if (term != TERM_HOPPING || !BaseType::hasSpinOrbit())  {
			edofs[0] = edofs[1] = 0;
			return;
		}

		SizeType spin1 = (dofs & 1);
		SizeType spin2 = (dofs & 2);
		spin2 >>= 1;

		// spin dependence of the hopping parameter (spin orbit)
		edofs[0] = spin1;
		edofs[1] = spin2;
	}

	SizeType terms() const { return 5; }
}; // class LinkProdExtendedSuperHubbard1Orb
} // namespace Dmrg
/*@}*/
#endif // LINKPROD_EXTENDED_SUPER_HUBBARD_1ORB_H
