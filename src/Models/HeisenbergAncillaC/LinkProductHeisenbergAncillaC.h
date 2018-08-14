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

/*! \file LinkProductHeisenbergAncillaC.h
 *
 *  LinkProduct for Heisenberg model
 *
 */
#ifndef DMRG_LINK_PROD_HEISENBERG_ANCILLAC_H
#define DMRG_LINK_PROD_HEISENBERG_ANCILLAC_H
#include "ProgramGlobals.h"
#include "LinkProductBase.h"

namespace Dmrg {
template<typename ModelHelperType, typename GeometryType>
class LinkProductHeisenbergAncillaC : public LinkProductBase<ModelHelperType, GeometryType> {

	typedef LinkProductBase<ModelHelperType, GeometryType> BaseType;
	typedef typename BaseType::AdditionalDataType AdditionalDataType;
	typedef typename BaseType::VectorSizeType VectorSizeType;

public:

	enum {TERM_SPLUSSMINUS, TERM_SZSZ, TERM_ANCILLA};

	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef typename ModelHelperType::RealType RealType;

	template<typename SomeInputType>
	LinkProductHeisenbergAncillaC(SomeInputType& io, bool hot)
	    : BaseType(io, "SplusiSminusj SziSzj Ancilla"), hot_(hot)
	{}

	void setLinkData(SizeType term,
	                 SizeType dofs,
	                 bool isSu2,
	                 ProgramGlobals::FermionOrBosonEnum& fermionOrBoson,
	                 std::pair<SizeType,SizeType>& ops,
	                 std::pair<char,char>&,
	                 SizeType& angularMomentum,
	                 RealType& angularFactor,
	                 SizeType& category,
	                 const AdditionalDataType&) const
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

	void valueModifier(SparseElementType& value,
	                   SizeType term,
	                   SizeType,
	                   bool isSu2,
	                   const AdditionalDataType&) const
	{
		if (term == 0) value *= 0.5;

		if (isSu2 && term != TERM_ANCILLA)
			value = -value;
	}

	SizeType dofs(SizeType term,const AdditionalDataType&) const
	{
		return (!hot_ || term == TERM_ANCILLA) ? 1 : 2;
	}

	// has only dependence on orbital
	void connectorDofs(VectorSizeType& edofs,
	                   SizeType term,
	                   SizeType dofs,
	                   const AdditionalDataType&) const
	{
		if (!hot_ || term == TERM_ANCILLA)
			edofs[0] = edofs[1] = 0;
		else
			edofs[0] = edofs[1] = dofs;
	}

	//! Splus Sminus and
	//! Sz Sz
	//! delta^\dagger delta
	SizeType terms() const { return 3; }

private:

	PairType operatorDofs(SizeType term,bool isSu2) const
	{
		if (term == TERM_SPLUSSMINUS || term == TERM_ANCILLA)
			return PairType(term,term);
		SizeType x = (isSu2) ? 0 : 1;
		return PairType(x,x);
	}

	PairType operatorDofsHot(SizeType term,SizeType dofs,bool) const
	{
		if (term == TERM_ANCILLA) return PairType(4,4);
		assert(term == TERM_SPLUSSMINUS || term == TERM_SZSZ);
		SizeType offset = (term == TERM_SPLUSSMINUS) ? 0 : 2;
		return PairType(dofs+offset,dofs+offset);
	}

	bool hot_;
}; // class LinkProductHeisenbergAncillaC
} // namespace Dmrg
/*@}*/
#endif

