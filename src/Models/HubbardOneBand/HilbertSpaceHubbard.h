/*
Copyright (c) 2009-2013, UT-Battelle, LLC
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

/*! \file HilbertSpaceHubbard.h
 *
 *  This class represents the Hilbert space for the Hubbard Model
 *  States are represented with binary numbers. Two bits per site
 *  Bits meaning:
 *  00 0 empty state
 *  01 1 up electron
 *  10 2 down electron
 *  11 3 doubly-occupied state
 * 
 *  Note: this is a static class
 * 
 */
#ifndef HILBERTSPACEHUBBARD_HEADER_H
#define HILBERTSPACEHUBBARD_HEADER_H


namespace Dmrg {

	//! A class to operate on quaterny numbers (base 4)
	template<typename Word>
	class HilbertSpaceHubbard {
	public:
		typedef Word HilbertState;	
		static int const  SPIN_UP=0;
		static int const  SPIN_DOWN=1;
		
		//! For state "a" set electron on site "j" to value "value"
		static void set(Word &a,int j,int value) 
		{
			Word mask; 
			switch (value) {
				case 0:
					mask = (1<<(2*j)) | (1<<(2*j+1));
					a &= (~mask);
					return;
				case 1:
					mask = (1<<(2*j));
					a |= mask;
					mask = (1<<(2*j+1));
					a &= (~mask);
					return;
				case 2:
					mask = (1<<(2*j+1));
					a |= mask;
					mask = (1<<(2*j));
					a &= (~mask);
					return;
				case 3:
					mask = (1<<(2*j)) | (1<<(2*j+1));
					a |= (mask);
					return;
				default:
					std::cerr<<"value="<<value<<"\n";
					throw PsimagLite::RuntimeError("HilbertSpaceHubbard: set: invalid value.\n");
					
						
			}
		}
		
		// Get electronic state on site "j" in binary number "a"
		static int get(Word const &a,int j)
		{
			Word mask = (1<<(2*j)) | (1<<(2*j+1));
			mask &= a;
			mask >>= (2*j);
			if (mask>3) {
				std::cerr<<"Error: mask="<<mask<<"\n";
			 	throw PsimagLite::RuntimeError("HilbertSpaceHubbard: get: invalid  mask.\n");
			}
			return mask;
			
		}
		
		// Destroy electron with internal dof  "sigma" on site "j" in binary number "a"
		static void destroy(Word &a,int j,int sigma)
		{
			Word mask;
			switch (sigma) {
				case 0:
					mask = (1<<(2*j));
					a &= (~mask);
					return;
				case 1:
					mask = (1<<(2*j+1));
					a &= (~mask);
					return;
				default:
					std::cerr<<"sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("HilbertSpaceHubbard: destroy: invalid value.\n");		
			}
		}
		
		// Create electron with internal dof  "sigma" on site "j" in binary number "a"
		static void create(Word &a,int j,int sigma)
		{
			Word mask;
			switch (sigma) {
				case 0:
					mask = (1<<(2*j));
					a |= mask;
					return;
				case 1:
					mask = (1<<(2*j+1));
					a |= mask;
					return;
				default:
					std::cerr<<"sigma="<<sigma<<"\n";
					throw PsimagLite::RuntimeError("HilbertSpaceHubbard: create: invalid value.\n");		
			}
		}
		
		// Is there an electron with internal dof  "sigma" on site "i" in binary number "ket"?
		static bool isNonZero(Word const &ket,int i,int sigma)
		{
			
			int tmp=get(ket,i);
			if ((tmp & 1) && sigma==0) return true;
			if ((tmp & 2) && sigma==1) return true;
			
			return false;
		}
		
		//! returns the number of electrons of internal dof "value" in binary number "data"
		static int getNofDigits(Word const &data,int value)
		{
			int ret=0;
			Word data2=data;
			int i=0;
			do {
				 if ( (data & (1<<(2*i+value))) ) ret++;
				 i++;
			} while (data2>>=1);
			
			return ret;
		}
		
		//! Number of electrons with dof sector between i and j excluding i and j in binary number "ket"
		//!  intended for when i<j
		static int calcNofElectrons(Word const &ket,int i,int j,int sector)
		{
			int ii=i+1;
			if (ii>=j) return 0;
			Word m=0;
			for (int k=2*ii;k<2*j;k++) m |= (1<<k);
			m = m & ket;
			return getNofDigits(m,sector);
		}
		
		
	}; // class HilbertSpaceHubbard
} // namespace Dmrg

/*@}*/	
#endif
