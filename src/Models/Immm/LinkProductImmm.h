// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009-2011, UT-Battelle, LLC
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

/*! \file LinkProductImmm.h
 *
 *  A class to represent product of operators that form a link or bond for this model
 *
 */
#ifndef LINK_PRODUCT_IMMM_H
#define LINK_PRODUCT_IMMM_H

namespace Dmrg {
	
	template<typename ModelHelperType>
	class LinkProductImmm {
			typedef typename ModelHelperType::SparseMatrixType SparseMatrixType;
			typedef typename SparseMatrixType::value_type SparseElementType;
			typedef std::pair<size_t,size_t> PairType;

		public:

			typedef typename ModelHelperType::RealType RealType;

			LinkProductImmm(size_t linSize)
			: linSize_(linSize),degreesOfFreedom_(linSize)
			{
				size_t cooperEach = 4;
				buildDofs(cooperEach);
			}
			
			size_t dOf(size_t i) const { return degreesOfFreedom_[i]; }

			void updateSites(size_t actualSite1,size_t actualSite2) const
			{
				actualSite1_=actualSite1;
				actualSite2_=actualSite2;
			}
			
			//! There are 4 different orbitals
			//! and 2 spins. Spin is diagonal so we end up with 8 possiblities
			//! a up a up, a up b up, b up a up, b up, b up,
			//! and similarly for spin down.
			size_t dofs(size_t term) const
			{ 
				//! Need to think this more, this isn't optimal
				//! because dofs(term) varies from site to site
				return 8;
			}

			// has only dependence on orbital
			PairType connectorDofs(size_t term,size_t dofs) const
			{
				//! Need to think this more
				//! This depends on the site
				size_t spin = dofs/4;
				size_t xtmp = (spin==0) ? 0 : 4;
				xtmp = dofs - xtmp;
				size_t orb1 = xtmp/2;
				size_t orb2 = (xtmp & 1);
				return PairType(orb1,orb2); // has only dependence on orbital
			}

			void setLinkData(size_t term,size_t dofs,bool isSu2,size_t& fermionOrBoson,PairType& ops,
     					std::pair<char,char>& mods,
					size_t& angularMomentum,
     					RealType& angularFactor,
					size_t& category) const
			{
				//!FIXME: Depends on site!!!!
				fermionOrBoson = ProgramGlobals::FERMION;
				size_t spin = getSpin(dofs);
				ops = operatorDofs(dofs);
				angularFactor = 1;
				if (spin==1) angularFactor = -1;
				angularMomentum = 1;
				category = spin;
			}
			
			void valueModifier(SparseElementType& value,size_t term,size_t dofs,bool isSu2) const
			{}

		private:

			// spin is diagonal
			std::pair<size_t,size_t> operatorDofs(size_t dofs) const
			{
				//!FIXME: Depends on site!!!!
				size_t spin = dofs/4;
				size_t xtmp = (spin==0) ? 0 : 4;
				xtmp = dofs - xtmp;
				size_t orb1 = xtmp/2;
				size_t orb2 = (xtmp & 1);
				size_t op1 = orb1 + spin*2;
				size_t op2 = orb2 + spin*2;
				return std::pair<size_t,size_t>(op1,op2);
			}
			
			size_t getSpin(size_t dofs) const
			{
				//!FIXME: Depends on site!!!!
				return dofs/4;
			}
			
			//! If there's only spin  at site i degreesOfFreedom_[i]=2
			//! If there's spin an norb orbitals then degreesOfFreedom_[i]=2*norb
			//! etc.
			void buildDofs(size_t copperEach)
			{
				size_t counter = 5;
				for (size_t i=0;i<degreesOfFreedom_.size();i++) {
					if (counter%copperEach==0) degreesOfFreedom_[i]=2;
					else degreesOfFreedom_[i]=4;
					counter++;
				}
			}

			size_t linSize_;
			std::vector<size_t> degreesOfFreedom_;
			mutable size_t actualSite1_;
			mutable size_t actualSite2_;
	}; // class LinkProductImmm
} // namespace Dmrg
/*@}*/
#endif // LINK_PRODUCT_IMMM_H
