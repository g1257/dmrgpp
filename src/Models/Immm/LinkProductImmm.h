/*
Copyright (c) 2009-2012, UT-Battelle, LLC
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

/*! \file LinkProductImmm.h
 *
 *  A class to represent product of operators that form a link or bond for this model
 *
 */
#ifndef LINK_PRODUCT_IMMM_H
#define LINK_PRODUCT_IMMM_H
#include <cassert>
#include "ProgramGlobals.h"
#include "LinkProductBase.h"

namespace Dmrg {

template<typename ModelHelperType>
class LinkProductImmm : public LinkProductBase<ModelHelperType> {

	typedef LinkProductBase<ModelHelperType> BaseType;
	typedef BaseType::AdditionalDataType AdditionalDataType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef std::pair<SizeType,SizeType> PairType;

	enum {HOPPING_TERM,W_TERM}; // W_TERM is a term of the form Upd n_i n_j

public:

	typedef typename ModelHelperType::RealType RealType;

	//! The term=0 is for hoppings:
	//! if sites are both TYPE_O then return 8, 2 orbitals, 2 spins,
	//! i.e. gamma sigma gamma' sigma = 2x2x2 = 8
	//! if sites are TYPE_O and TYPE_C the returns 4, i.e. gamma sigma d sigma =  2x2
	//! The term=1 is for Upd
	SizeType dofs(SizeType term,const AdditionalDataType& additionalData)
	{
		if (term==W_TERM) {
			// No Upd allowed between Oxygens
			if (additionalData.type1==additionalData.type2)
				return 0;
			return 1; // Upd n_o n_Cu
		}
		SizeType type1 = additionalData.type1;
		SizeType type2 = additionalData.type2;
		//! both cannot be TYPE_C
		assert(type1!=type2 || type1!=additionalData.TYPE_C);
		return (type1==type2) ? 8 : 4;
	}

	// has only dependence on orbital
	void connectorDofs(VectorSizeType& edofs,
	                   SizeType term,
	                   SizeType dofs,
	                   const AdditionalDataType& additionalData)
	{
		if (term==W_TERM) {
			edofs[0] = edofs[1] = 0;
			return;
		}; // Upd

		SizeType type1 = additionalData.type1;
		SizeType type2 = additionalData.type2;
		//! both cannot be TYPE_C
		assert(type1!=type2 || type1!=additionalData.TYPE_C);
		//! both TYPE_O:
		if (type1==type2) {
			SizeType spin = dofs/4;
			SizeType xtmp = (spin==0) ? 0 : 4;
			xtmp = dofs - xtmp;
			SizeType orb1 = xtmp/2;
			SizeType orb2 = (xtmp & 1);
			edofs[0] = orb1;
			edofs[1] = orb2;
			return; // has only dependence on orbital
		}

		//! TYPE_C and TYPE_O:
		SizeType spin = dofs/2;
		SizeType xtmp = (spin==0) ? 0 : 2;
		if (type1==additionalData.TYPE_C) {
			edofs[0] = 0;
			edofs[1] = dofs - xtmp;
		} else {
			edofs[0] = dofs - xtmp;
			edofs[1] = 0;
		}
	}

	void setLinkData(SizeType term,
	                 SizeType dofs,
	                 bool,
	                 ProgramGlobals::FermionOrBosonEnum& fermionOrBoson,
	                 PairType& ops,
	                 std::pair<char,char>&,
	                 SizeType& angularMomentum,
	                 RealType& angularFactor,
	                 SizeType& category,
	                 const AdditionalDataType& additionalData)
	{
		if (term==W_TERM) {
			fermionOrBoson = ProgramGlobals::BOSON;
			assert(dofs==0);
			if (additionalData.type1==additionalData.TYPE_C) {
				ops.first = 2;
				ops.second = 4;
			} else {
				ops.first = 4;
				ops.second = 2;
			}
			angularFactor = 1;
			angularMomentum = 1;
			category = 0;
			return;
		}

		fermionOrBoson = ProgramGlobals::FERMION;
		SizeType spin = getSpin(dofs,additionalData);
		ops = operatorDofs(dofs,additionalData);
		angularFactor = 1;
		if (spin==1) angularFactor = -1;
		angularMomentum = 1;
		category = spin;
	}

	void valueModifier(SparseElementType& value,
	                   SizeType term,
	                   SizeType,
	                   bool,
	                   const AdditionalDataType&)
	{
		if (term==W_TERM) {
			value *= 0.5;
			return;
		}

		value *= (-1.0);
	}

	SizeType terms() { return 2; }

private:

	// spin is diagonal
	std::pair<SizeType,SizeType> operatorDofs(SizeType dofs,
	                                          const AdditionalDataType& additionalData)
	{
		SizeType type1 = additionalData.type1;
		SizeType type2 = additionalData.type2;
		//! both cannot be TYPE_C
		assert(type1!=type2 || type1!=additionalData.TYPE_C);
		//! both TYPE_O:
		if (type1==type2) {
			SizeType spin = dofs/4;
			SizeType xtmp = (spin==0) ? 0 : 4;
			xtmp = dofs - xtmp;
			SizeType orb1 = xtmp/2;
			SizeType orb2 = (xtmp & 1);
			SizeType op1 = orb1 + spin*2;
			SizeType op2 = orb2 + spin*2;
			assert(op1<4 && op2<4);
			return std::pair<SizeType,SizeType>(op1,op2);
		}

		//! TYPE_C and TYPE_O:
		assert(dofs<4);
		SizeType spin = dofs/2;
		SizeType xtmp = (spin==0) ? 0 : 2;
		xtmp = dofs - xtmp;
		SizeType op1 = spin;
		SizeType op2 = xtmp + spin*2;
		assert(op1<2 && op2<4);
		return (type1==additionalData.TYPE_C) ? PairType(op1,op2) : PairType(op2,op1);
	}

	SizeType getSpin(SizeType dofs,const AdditionalDataType& additionalData)
	{
		SizeType type1 = additionalData.type1;
		SizeType type2 = additionalData.type2;
		//! both cannot be TYPE_C
		assert(type1!=type2 || type1!=additionalData.TYPE_C);

		return (type1==type2) ? dofs/4 : dofs/2;
	}
}; // class LinkProductImmm
} // namespace Dmrg
/*@}*/
#endif // LINK_PRODUCT_IMMM_H

