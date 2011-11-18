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

/*! \file HilbertSpaceImmm.h
 *
 *  This class represents the Hilbert space for the Immm Model
 *  States are represented with binary numbers. N(site) bits per site
 *  Bits meaning:
 *  0.....0  empty state
 *  0..1..0  a 1 at location x means state "x"
 *  0..010..010..0, a 1 at location x and a 1 at location y means 2 electrons on the site with states x and y respectively 
 *  ...
 *  1111...111111  all ones means N electrons each with a different state
 * 
 *  Note: this is a class
 *  Note: Length of state depends on site
 * 
 */
#ifndef HILBERTSPACE_IMMM_H
#define HILBERTSPACE_IMMM_H
#include <iostream>
#include <vector>

namespace Dmrg {

	//! A class to operate on n-ary numbers (base n)
	template<typename Word,typename LinkProductType>
	class HilbertSpaceImmm {

	public:

		typedef Word HilbertState;

		static const size_t NUMBER_OF_SPINS = 2;
		enum {SPIN_UP=0,SPIN_DOWN=1};

		HilbertSpaceImmm(const LinkProductType& linkProduct)
		: linkProduct_(linkProduct)
		{}

		size_t dOf(size_t i) const { return linkProduct_.dOf(i); }
		
		// Get electronic state on site "j" in binary number "a"
		Word get(Word const &a,size_t j) const
		{
			size_t k=degreesOfFreedomUpTo(j);
			size_t ones = (1<<(dOf(j)))-1;
			Word mask=(ones<<k);

			mask &= a;
			mask >>= k;
			return mask;
			
		}

		// Create electron with internal dof  "sigma" on site "j" in binary number "a"
		void create(Word &a,size_t j,size_t sigma) const
		{
			size_t k=degreesOfFreedomUpTo(j);
			Word mask=(1<<(k+sigma));
			a |= mask;
		}
		
		// Is there an electron with internal dof  "sigma" on site "i" in binary number "ket"?
		bool isNonZero(Word const &ket,size_t i,size_t sigma) const
		{
			
			Word tmp=get(ket,i);
			//std::cerr<<"isNonZero, ket="<<ket<<" tmp="<<tmp<<"\n";
			if (tmp & (1<<sigma)) return true;
			
			return false;
		}
		
		//! returns the number of electrons of internal dof "value" in binary number "data"
		int getNofDigits(const Word& data,size_t value) const
		{
			int ret=0;
			Word data2=data;
			size_t i=0;
			size_t dof = 0;

			do {
				size_t k=degreesOfFreedomUpTo(i);
				dof = dOf(i);
				if ( (data & (1<<(k+value))) ) ret++;
				i++;
			} while (data2>>=dof);
			
			return ret;
		}
		
		//! Number of electrons with spin spin (sums over bands)
		int electronsWithGivenSpin(Word const &data,size_t site,size_t spin) const
		{
			
			size_t norb = dOf(site)/NUMBER_OF_SPINS;
			size_t beginX=spin*norb;
			size_t endX=beginX + norb;
			size_t sum=0;

			for (size_t x=beginX;x<endX;x++) sum += getNofDigits(data,x);
			
			return sum;	
			
		}

		//! Number of electrons at given site (sum over all bands)
		int electronsAtGivenSite(const Word& data,size_t site) const
		{
			size_t sum=0;
			
			for (size_t sector=0;sector<dOf(site);sector++)
				sum += calcNofElectrons(data,0,sector);
			
			return sum;	
			
		}

		//! Number of electrons with dof sector between i and j excluding i and j in binary number "ket"
		//!  intended for when i<j
		 int calcNofElectrons(const Word& ket,size_t i,size_t j,size_t sector) const
		{
			size_t ii=i+1;
			if (ii>=j) return 0;
			Word m=0;
			for (size_t site = ii;site<j;site++) {
				size_t k = degreesOfFreedomUpTo(site);
				size_t dof = dOf(site);
				for (size_t sigma=0;sigma<dof;sigma++)
					m |= (1<<(k+sigma));
			}
			m &= ket;
			return getNofDigits(m,sector);
		} 

		//! Number of electrons with dof sector on site i in binary number "ket"
		int calcNofElectrons(Word const &ket,size_t i,size_t sector) const
		{
			Word m=0;
			
			size_t k = degreesOfFreedomUpTo(i);
			size_t dof = dOf(i);
			for (size_t sigma=0;sigma<dof;sigma++) m |= (1<<(k+sigma));
			
			m &= ket;
			return getNofDigits(m,sector);
		}

	private:

		size_t degreesOfFreedomUpTo(size_t j) const
		{
			size_t k=0;
			for (size_t kk=0;kk<j;kk++) k += dOf(kk);
			return k;
		}

		const LinkProductType& linkProduct_;
	}; // class HilbertSpaceImmm
} // namespace Dmrg

/*@}*/	
#endif // HILBERTSPACE_IMMM_H
