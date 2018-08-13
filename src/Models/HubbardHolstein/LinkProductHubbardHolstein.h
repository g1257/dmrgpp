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

/*! \file LinkProductHubbardHolstein.h
 *
 *  A class to represent product of operators that form a link or
 *  bond for this model
 *
 */
#ifndef DMRG_LINKPROD_HUBBARD_HOLSTEIN_H
#define DMRG_LINKPROD_HUBBARD_HOLSTEIN_H
#include "ProgramGlobals.h"
#include "LinkProductBase.h"

namespace Dmrg {

template<typename ModelHelperType, typename GeometryType>
class LinkProductHubbardHolstein : public LinkProductBase<ModelHelperType, GeometryType> {

	typedef LinkProductBase<ModelHelperType, GeometryType> BaseType;
	typedef typename GeometryType::AdditionalDataType AdditionalDataType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef typename SparseMatrixType::value_type SparseElementType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef typename ModelHelperType::BasisType BasisType;
	typedef typename ModelHelperType::OperatorType OperatorType;

	enum {TERM_HOPPINGF, TERM_HOPPINGP, TERM_HOPPINGSSH};

public:

	typedef std::pair<char,char> PairCharType;
	typedef typename ModelHelperType::RealType RealType;

	template<typename SomeInputType>
	LinkProductHubbardHolstein(SomeInputType& io, bool isSsh)
	    : BaseType(io, (isSsh) ? "HoppingFermionic HoppingBosonic" :
	                             "HoppingFermionic HoppingBosonic HoppingSSH"),
	      isSsh_(isSsh),
	      numberphonons_(0)
	{
		io.readline(numberphonons_,"NumberPhonons=");
	}

	SizeType dofs(SizeType term, const AdditionalDataType&) const
	{
		if (term == TERM_HOPPINGF) return 2;
		if (term == TERM_HOPPINGP) return 1;
		assert(term == TERM_HOPPINGSSH);
		return 8;
	}

	void setLinkData(SizeType term,
	                 SizeType dof,
	                 bool,
	                 ProgramGlobals::FermionOrBosonEnum& fermionOrBoson,
	                 PairType& ops,
	                 PairCharType& mods,
	                 SizeType& angularMomentum,
	                 RealType& angularFactor,
	                 SizeType& category,
	                 const AdditionalDataType&) const
	{
		mods = PairCharType('C', 'N');

		if (term==TERM_HOPPINGF) {
			assert(dof == 0 || dof == 1);
			fermionOrBoson = ProgramGlobals::FERMION;
			ops = PairType(dof, dof);
			angularFactor = 1;
			if (dof == 1) angularFactor = -1;
			angularMomentum = 1;
			category = dof;
		}

		if (term==TERM_HOPPINGP) {
			assert(dof == 0);
			fermionOrBoson = ProgramGlobals::BOSON;
			SizeType offset1 = 2;
			ops = PairType(dof + offset1, dof + offset1);
		}

		if (term==TERM_HOPPINGSSH) {
			mods = PairCharType('C', 'N');
			fermionOrBoson = ProgramGlobals::FERMION;
			SizeType offset2 = 3;
			assert(dof >= 0 && dof < 8);
			switch (dof) {
			case 0: // old 0
				ops = PairType(0, offset2);
				break;
			case 1: // old 1
				ops = PairType(offset2, 0);
				break;
			case 2: // old 2
				ops = PairType(1, offset2 + 1);
				break;
			case 3: // old 3
				ops = PairType(offset2 + 1, 1);
				break;
			}

			angularFactor = 1;
			if (dof == 1) angularFactor = -1;
			angularMomentum = 1;
			category = dof;
		}
	}

	SizeType terms() const { return (isSsh_) ? 3 : 2; }

private:

	bool isSsh_;
	SizeType numberphonons_;
}; // class LinkProductHubbardHolstein
} // namespace Dmrg
/*@}*/
#endif //DMRG_LINKPROD_HUBBARD_HOLSTEIN_H

