/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file LinkProdSuperHubbardExtended.h
 *
 *  FIXME
 *
 */
#ifndef LINKPROD_SUPER_HUBBARD_EXTENDED_H
#define LINKPROD_SUPER_HUBBARD_EXTENDED_H

#include "../Models/Heisenberg/LinkProductHeisenberg.h"
#include "../Models/ExtendedHubbard1Orb/LinkProdExtendedHubbard1Orb.h"
#include "LinkProductBase.h"

namespace Dmrg {

template<typename ModelHelperType, typename GeometryType>
class LinkProdSuperHubbardExtended : public LinkProductBase<ModelHelperType, GeometryType> {

	typedef LinkProductBase<ModelHelperType, GeometryType> BaseType;
	typedef typename BaseType::AdditionalDataType AdditionalDataType;
	typedef typename BaseType::VectorSizeType VectorSizeType;
	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef std::pair<SizeType,SizeType> PairType;
	typedef LinkProdExtendedHubbard1Orb<ModelHelperType,
	GeometryType> LinkProdExtendedHubbard1OrbType;
	typedef LinkProductHeisenberg<ModelHelperType, GeometryType> LinkProductHeisenbergType;

	enum {TERM_HOPPING=0 ,TERM_NINJ=1, TERM_SUPER = 2};

	static const bool ANISOTROPIC_IS_FALSE = false;

public:

	typedef typename ModelHelperType::RealType RealType;
	typedef typename SparseMatrixType::value_type SparseElementType;

	template<typename SomeInputType>
	LinkProdSuperHubbardExtended(SomeInputType& io)
	    : BaseType(io, "Hopping NiNj Super"),
	      lHeis_(io, ANISOTROPIC_IS_FALSE),
	      lExtHubb_(io)
	{}

	void setLinkData(SizeType term,
	                 SizeType dofs,
	                 bool isSu2,
	                 ProgramGlobals::FermionOrBosonEnum& fermionOrBoson,
	                 PairType& ops,
	                 std::pair<char,char>& mods,
	                 SizeType& angularMomentum,
	                 RealType& angularFactor,
	                 SizeType& category,
	                 const AdditionalDataType& additional) const
	{
		if (term == TERM_HOPPING || term == TERM_NINJ) {
			return lExtHubb_.setLinkData(term,
			                             dofs,
			                             isSu2,
			                             fermionOrBoson,
			                             ops,
			                             mods,
			                             angularMomentum,
			                             angularFactor,
			                             category,
			                             additional);
		}

		assert(term == TERM_SUPER);

		lHeis_.setLinkData(0,
		                   dofs,
		                   isSu2,
		                   fermionOrBoson,
		                   ops,
		                   mods,
		                   angularMomentum,
		                   angularFactor,
		                   category,
		                   additional);
		ops.first += 3;
		ops.second += 3;
	}

	void valueModifier(SparseElementType& value,
	                   SizeType term,
	                   SizeType dofs,
	                   bool isSu2,
	                   const AdditionalDataType& additional) const
	{
		if (term == TERM_HOPPING || term == TERM_NINJ)
			return lExtHubb_.valueModifier(value,
			                               term,
			                               dofs,
			                               isSu2,
			                               additional);

		return lHeis_.valueModifier(value,
		                            0,
		                            dofs,
		                            isSu2,
		                            additional);
	}

	SizeType dofs(SizeType term,const AdditionalDataType& additional) const
	{
		if (term == TERM_HOPPING || term == TERM_NINJ)
			return lExtHubb_.dofs(term,
			                      additional);

		return lHeis_.dofs(0,additional);
	}

	// has only dependence on orbital
	void connectorDofs(VectorSizeType& edofs,
	                   SizeType term,
	                   SizeType dofs,
	                   const AdditionalDataType& additional) const
	{
		if (term == TERM_HOPPING || term == TERM_NINJ)
			return lExtHubb_.connectorDofs(edofs,
			                               term,
			                               dofs,
			                               additional);

		return lHeis_.connectorDofs(edofs,
		                            0,
		                            dofs,
		                            additional);
	}

	SizeType terms() const { return 3; }

private:

	LinkProductHeisenbergType lHeis_;
	LinkProdExtendedHubbard1OrbType lExtHubb_;
}; // class LinkProdSuperHubbardExtended
} // namespace Dmrg
/*@}*/
#endif // LINKPROD_SUPER_HUBBARD_EXTENDED_H
