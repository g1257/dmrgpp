/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file LinkProdSuperHubbardExtended.h
 *
 *  FIXME
 *
 */
#ifndef LINKPROD_SUPER_HUBBARD_EXTENDED_H
#define LINKPROD_SUPER_HUBBARD_EXTENDED_H

#include "../Models/Heisenberg/LinkProductHeisenberg.h"
#include "../Models/ExtendedHubbard1Orb/LinkProdExtendedHubbard1Orb.h"

namespace Dmrg {

template<typename ModelHelperType>
class LinkProdSuperHubbardExtended {

	typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
	typedef std::pair<SizeType,SizeType> PairType;
	enum {TERM_HOPPING=0 ,TERM_NINJ=1, TERM_SUPER = 2};

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
	                        const SomeStructType& additional)
	{
		if (term == TERM_HOPPING || term == TERM_NINJ) {
			return LinkProdExtendedHubbard1Orb<ModelHelperType>::setLinkData(term,
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

		LinkProductHeisenberg<ModelHelperType>::setLinkData(0,
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

	template<typename SomeStructType>
	static void valueModifier(SparseElementType& value,
	                          SizeType term,
	                          SizeType dofs,
	                          bool isSu2,
	                          const SomeStructType& additional)
	{
		if (term == TERM_HOPPING || term == TERM_NINJ)
			return LinkProdExtendedHubbard1Orb<ModelHelperType>::valueModifier(value,
			                                                                   term,
			                                                                   dofs,
			                                                                   isSu2,
			                                                                   additional);

		return LinkProductHeisenberg<ModelHelperType>::valueModifier(value,
		                                                             0,
		                                                             dofs,
		                                                             isSu2,
		                                                             additional);
	}

	template<typename SomeStructType>
	static SizeType dofs(SizeType term,const SomeStructType& additional)
	{
		if (term == TERM_HOPPING || term == TERM_NINJ)
			return LinkProdExtendedHubbard1Orb<ModelHelperType>::dofs(term,additional);

		return LinkProductHeisenberg<ModelHelperType>::dofs(0,additional);
	}

	template<typename SomeStructType>
	static PairType connectorDofs(SizeType term,
	                              SizeType dofs,
	                              const SomeStructType& additional)
	{
		if (term == TERM_HOPPING || term == TERM_NINJ)
			return LinkProdExtendedHubbard1Orb<ModelHelperType>::connectorDofs(term,
			                                                                   dofs,
			                                                                   additional);

		return LinkProductHeisenberg<ModelHelperType>::connectorDofs(0,
		                                                             dofs,
		                                                             additional);
	}

	static SizeType terms() { return 3; }

}; // class LinkProdSuperHubbardExtended
} // namespace Dmrg
/*@}*/
#endif // LINKPROD_SUPER_HUBBARD_EXTENDED_H
