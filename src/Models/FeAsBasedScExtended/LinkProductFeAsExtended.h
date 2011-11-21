// BEGIN LICENSE BLOCK
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
// END LICENSE BLOCK
/** \ingroup DMRG */
/*@{*/

/*! \file LinkProductFeAsExtended.h
 *
 *  A class to represent product of operators that form a link or
 *  bond for this model
 *
 */
#ifndef LINK_PRODUCT_EX_H
#define LINK_PRODUCT_EX_H
#include "LinkProductHeisenbergSpinOneHalf.h"
#include "LinkProductFeAs.h"
#include "OperatorsFeAs.h"

namespace Dmrg {
	
	template<typename ModelHelperType>
	class LinkProductFeAsExtended {
			typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef std::pair<size_t,size_t> PairType;
			typedef LinkProductFeAs<ModelHelperType> LinkProductFeAsType;
			typedef LinkProductHeisenbergSpinOneHalf<ModelHelperType>
				LinkProductHeisenbergType;
			typedef typename ModelHelperType::BasisType BasisType;
			typedef typename ModelHelperType::OperatorType OperatorType;
			typedef OperatorsFeAs<OperatorType,BasisType> OperatorsFeAsType;

			enum {TERM_HOPPING, TERM_J};
			
		public:
			typedef typename ModelHelperType::RealType RealType;
			
			template<typename SomeStructType>
			static size_t dofs(size_t term,const SomeStructType& additional)
			{
				return (term==TERM_J) ?
						LinkProductHeisenbergType::dofs(term,additional) :
						LinkProductFeAsType::dofs(term,additional);
			}
			
			// has only dependence on orbital
			template<typename SomeStructType>
			static PairType connectorDofs(size_t term,size_t dofs,const SomeStructType& additional)
			{
				if (term==TERM_HOPPING)
					return LinkProductFeAsType::connectorDofs(term,dofs,additional);

				return LinkProductHeisenbergType::connectorDofs(term,dofs,additional);
			}

			template<typename SomeStructType>
			static void setLinkData(
					size_t term,
					size_t dofs,
     					bool isSu2,
					size_t& fermionOrBoson,
					PairType& ops,
     					std::pair<char,char>& mods,
					size_t& angularMomentum,
     					RealType& angularFactor,
					size_t& category,const SomeStructType& additional)
			{
				if (term==TERM_HOPPING)
					return LinkProductFeAsType::setLinkData(
						term,dofs,isSu2,fermionOrBoson,ops,mods,
						angularMomentum,angularFactor,category,additional);

				LinkProductHeisenbergType::setLinkData(
						term,dofs,isSu2,fermionOrBoson,ops,mods,
						angularMomentum,angularFactor,category,additional);
				size_t offset1 = OperatorsFeAsType::DEGREES_OF_FREEDOM;
				ops.first += offset1;
				ops.second += offset1;
			}

			template<typename SomeStructType>
			static void valueModifier(
					SparseElementType& value,
					size_t term,
					size_t dofs,bool isSu2,const SomeStructType& additional)
			{
				if (term==TERM_HOPPING) return
						LinkProductFeAsType::valueModifier(
							value,term,dofs,isSu2,additional);
				LinkProductHeisenbergType::valueModifier(
						value,term,dofs,isSu2,additional);

			}
	}; // class LinkProductFeAsExtended
} // namespace Dmrg
/*@}*/
#endif //LINK_PRODUCT_EX_H
