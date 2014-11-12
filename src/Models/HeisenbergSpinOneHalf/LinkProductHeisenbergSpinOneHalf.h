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

/*! \file LinkProductHeisenbergSpinOneHalf.h
 *
 *  LinkProduct for HeisenbergSpinOneHalf model
 *
 */
#ifndef LINK_PRODUCT_HEIS_ONE_HALF_H
#define LINK_PRODUCT_HEIS_ONE_HALF_H

namespace Dmrg {
	template<typename ModelHelperType>
	class LinkProductHeisenbergSpinOneHalf {

		public:
			typedef std::pair<SizeType,SizeType> PairType;
			typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef typename ModelHelperType::RealType RealType;

			template<typename SomeStructType>
			static void setLinkData(
					SizeType,
					SizeType dofs,
     					bool isSu2,
					SizeType& fermionOrBoson,
					std::pair<SizeType,SizeType>& ops,
					std::pair<char,char>& mods,
					SizeType& angularMomentum,
     					RealType& angularFactor,
					SizeType& category,
					const SomeStructType&)
			{
				fermionOrBoson = ProgramGlobals::BOSON;
				ops = operatorDofs(dofs,isSu2);
				angularMomentum = 2;
				char tmp = mods.first;
				switch (dofs) {
					case 0:
						angularFactor = -1;
						category = 2;
						break;
					case 1:
						angularFactor = -1;
						category = 0;
						mods.first = mods.second;
						mods.second = tmp;
						break;
					case 2:
						angularFactor = 0.5;
						category = 1;
						break;
				}
			}

			template<typename SomeStructType>
			static void valueModifier(SparseElementType& value,SizeType,SizeType dofs,bool isSu2,const SomeStructType&)
			{
				if (isSu2) value = -value;
				if (dofs<2) value *= 0.5;
				value *= 0.5;
			}

			//! For TERM_J there are 3 terms:
			//! Splus Sminus and
			//! Sminus Splus and
			//! Sz Sz
			template<typename SomeStructType>
			static SizeType dofs(SizeType,const SomeStructType&) { return 3; }

			template<typename SomeStructType>
			static PairType connectorDofs(SizeType,SizeType,const SomeStructType&)
			{
				return PairType(0,0); // no orbital and no anisotropy
			}

			static SizeType terms() { return 1; }

		private:

			static PairType operatorDofs(SizeType dofs,bool isSu2)
			{
				if (dofs<2) return PairType(0,0);
				SizeType x = (isSu2) ? 0 : 1;
				return PairType(x,x);
			}
	}; // class LinkProductHeisenbergSpinOneHalf
} // namespace Dmrg
/*@}*/
#endif
