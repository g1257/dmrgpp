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

/*! \file HilbertSpaceFeAs.h
 *
 *  This class represents the Hilbert space for the FeAs Model
 *  States are represented with binary numbers. N bits per site
 *  Bits meaning:
 *  0.....0  empty state
 *  0..1..0  a 1 at location x means state "x"
 *  0..010..010..0, a 1 at location x and a 1 at location y means 2 electrons on the site with states x and y respectively 
 *  ...
 *  1111...111111  all ones means N electrons each with a different state
 * 
 *  Note: this is a static class
 * 
 */
#ifndef HILBERTSPACEFEAS_HEADER_H
#define HILBERTSPACEFEAS_HEADER_H


namespace Dmrg {

	//! A class to operate on n-ary numbers (base n)
	template<typename Word>
	class HilbertSpaceFeAs {
	public:
		typedef Word HilbertState;	
		static size_t const NUMBER_OF_ORBITALS =2;
		static size_t const NUMBER_OF_STATES = 4; //2*NUMBER_OF_ORBITALS;
		enum {SPIN_UP=0,SPIN_DOWN=1};
		
		//! For state "a" set electron on site "j" to value "value"
		static void set(Word &a,size_t j,size_t value) 
		{
			
			size_t k=NUMBER_OF_STATES*j;
			size_t ones = 2*NUMBER_OF_STATES-1;
			Word mask=(ones<<k);
			Word b =(a & (~mask));
			
			mask=(value<<k);
			
			a= (b | mask);
		}

		
		// Get electronic state on site "j" in binary number "a"
		static Word get(Word const &a,size_t j)
		{
			size_t k=NUMBER_OF_STATES*j;
			size_t ones = (1<<(NUMBER_OF_STATES))-1;
			Word mask=(ones<<k);
			
			mask &= a;
			mask >>= k;
			return mask;
			
		}
		
		
		// Create electron with internal dof  "sigma" on site "j" in binary number "a"
		static void create(Word &a,size_t j,size_t sigma)
		{
			size_t k=NUMBER_OF_STATES*j;
			Word mask=(1<<(k+sigma));
			a |= mask;
		}
		
		// Is there an electron with internal dof  "sigma" on site "i" in binary number "ket"?
		static bool isNonZero(Word const &ket,size_t i,size_t sigma)
		{
			
			Word tmp=get(ket,i);
			//std::cerr<<"isNonZero, ket="<<ket<<" tmp="<<tmp<<"\n";
			if (tmp & (1<<sigma)) return true;
			
			return false;
		}
		
		//! returns the number of electrons of internal dof "value" in binary number "data"
		static int getNofDigits(Word const &data,size_t value)
		{
			int ret=0;
			Word data2=data;
			size_t i=0;
			do {
				 if ( (data & (1<<(NUMBER_OF_STATES*i+value))) ) ret++;
				 i++;
			} while (data2>>=NUMBER_OF_STATES);
			
			return ret;
		}
		
		//! Number of electrons with spin spin (sums over bands)
		static int electronsWithGivenSpin(Word const &data,size_t spin)
		{
			size_t sum=0;
			size_t beginX=spin*NUMBER_OF_ORBITALS;
			size_t endX=beginX + NUMBER_OF_ORBITALS;
			
			for (size_t x=beginX;x<endX;x++)
				sum += getNofDigits(data,x);
			
			return sum;	
			
		}
		
		//! Number of electrons at given site (sum over all bands)
		static int electronsAtGivenSite(Word const &data,size_t site)
		{
			size_t sum=0;
			
			for (size_t sector=0;sector<NUMBER_OF_STATES;sector++)
				sum += calcNofElectrons(data,site,sector);
			
			return sum;	
			
		}
		//! Number of electrons with dof sector between i and j excluding i and j in binary number "ket"
		//!  intended for when i<j
		 static int calcNofElectrons(Word const &ket,size_t i,size_t j,size_t sector)
		{
			size_t ii=i+1;
			if (ii>=j) return 0;
			Word m=0;
			for (size_t k=NUMBER_OF_STATES*ii;k<NUMBER_OF_STATES*j;k++) m |= (1<<k);
			m = m & ket;
			return getNofDigits(m,sector);
		} 
		
		//! Number of electrons with dof sector on site i in binary number "ket"
		 static int calcNofElectrons(Word const &ket,size_t i,size_t sector)
		{
			Word m=0;
			for (size_t k=NUMBER_OF_STATES*i;k<NUMBER_OF_STATES*(i+1);k++) m |= (1<<k);
			m = m & ket;
			return getNofDigits(m,sector);
		} 
		
	}; // class HilbertSpaceFeAs
} // namespace Dmrg

/*@}*/	
#endif
