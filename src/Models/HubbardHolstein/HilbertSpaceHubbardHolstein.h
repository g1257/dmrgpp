/*
Copyright (c) 2009-2014, UT-Battelle, LLC
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

/*! \file HilbertSpaceHubbardHolstein.h
 *
 *  This class represents the Hilbert space for the Hubbard Holstein Model
 *  States are represented with binary numbers. N bits per site (2 for the fermions, N-2 for the phonons)
 *  Bits meaning:
 *  0....00  empty state
 *  0..1.00  a 1 at location x>2 means state with N-2 maxphonons and zero fermions
 *  1111111  all ones means 2 electrons
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACE_HUBBARD_HOLSTEIN_HEADER_H
#define HILBERTSPACE_HUBBARD_HOLSTEIN_HEADER_H

#include "Utils.h"

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template<typename Word>
class HilbertSpaceHubbardHolstein {

	static SizeType numberBitphonons_;

public:

	typedef Word HilbertState;

	enum {SPIN_UP=0,SPIN_DOWN=1};

	static void setBitPhonons(SizeType numberphonons)
	{
		numberBitphonons_ =utils::bitSizeOfInteger(numberphonons);
	}

	// Get full electronic and phononic state on site "j" in binary number "a"
	static Word get(Word const &a,SizeType j)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType k=dofs*j;
		SizeType ones = (1<<(dofs))-1;
		Word mask=(ones<<k);

		mask &= a;
		mask >>= k;
		return mask;

	}

	// Get electronic state on site "j" in binary number "a"
	static Word getF(Word const &a,SizeType j)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType k=dofs*j;
		Word mask=(3<<k);

		mask &= a;
		mask >>= k;
		return mask;

	}

	// Get phononic state on site "j" in binary number "a"
	static Word getP(Word const &a,SizeType j)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType k=dofs*j;
		SizeType onesP = (~3);
		Word mask=(onesP<<k);

		mask &= a;
		mask >>= k;
		mask >>= 2; // Shitf by the Fermionic State
		return mask;

	}
	// Create electron with internal dof "sigma" on site "j" in binary number "a"
	static void createF(Word &a,SizeType j,SizeType sigma)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType k=dofs*j;
		Word mask=(1<<(k+sigma));
		a |= mask;
	}
	// Create phonon on site "j"
	static void createP(Word &a,SizeType j)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType k=dofs*j;
		SizeType nphonons = 0;
		if (isNonZeroP(a,j))
			nphonons = int(getP(a,j));

		Word stateP = 1+nphonons;
		Word maskP = (stateP<<2);
		Word maskF = getF(a,j);
		Word maskU = maskP | maskF;
		a = (maskU<<k);
	}

	// Destroy electron with internal dof "sigma" on site "j" in binary number "a"
	static void destroyF(Word &a,SizeType j,SizeType sigma)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType k=dofs*j;
		Word mask=(1<<(k+sigma));
		a &= (~mask);
	}

	// Is there a phonon on site "i" in binary number "ket"?
	static bool isNonZeroP(Word const &ket,SizeType i)
	{

		Word tmp=getP(ket,i);
		SizeType counter = 0;
		while (tmp) {
			if (tmp & 1) counter++;
			tmp >>=1;
		}
		if (counter) return true;

		return false;
	}

	// Is there an electron with internal dof "sigma" on site "i" in binary number "ket"?
	static bool isNonZeroF(Word const &ket,SizeType i,SizeType sigma)
	{

		Word tmp=getF(ket,i);
		if (tmp & (1<<sigma)) return true;

		return false;
	}

	//! returns the number of electrons of internal dof "sigma" in binary number "data"
	static int getNofDigits(Word const &data,SizeType sigma)
	{
		SizeType dofs = 2+numberBitphonons_;
		int ret=0;
		Word data2=data;
		SizeType i=0;
		do {
			if ( (data & (1<<(dofs*i+sigma))) ) ret++;
			i++;
		} while (data2>>=dofs);

		return ret;
	}

	//! Number of electrons with spin spin (sums over sites)
	static int electronsWithGivenSpin(Word const &data,SizeType spin)
	{
		SizeType sum=0;
		Word data2 = data;
		SizeType shift = 2+numberBitphonons_;
		while (data2 > 0) {
			if (data2 & (spin+1)) sum++;
			data2 >>= shift;
		}

		return sum;

	}

	//! Number of electrons in binary number "data" (sum over all spins)
	static SizeType electrons(const Word& data)
	{
		SizeType sum=0;
		SizeType dofs = 2;
		for (SizeType spin=0;spin<dofs;spin++)
			sum += calcNofElectrons(data,spin);

		return sum;

	}
	//! Number of electrons with spin spin between i and j excluding
	//! i and j in binary number "ket"
	//!  intended for when i<j

	static int calcNofElectrons(Word const &ket,SizeType i,SizeType j,SizeType spin)
	{
		SizeType dofs = 2+numberBitphonons_;
		SizeType ii=i+1;
		if (ii>=j) return 0;
		Word m=0;
		for (SizeType k=dofs*ii;k<dofs*j;k++) m |= (1<<k);
		m = m & ket;
		return getNofDigits(m,spin);
	}

	//! Number of electrons with spin spin in binary number "ket"
	static int calcNofElectrons(Word const &ket,SizeType spin)
	{
		SizeType dofs = 2+numberBitphonons_;
		Word ket2 = ket;
		SizeType digit = spin+1;
		int sum = 0;
		while (ket2 > 0) {
			if (ket2 & digit) sum++;
			ket2 >>= dofs;
		}

		return sum;
	}

}; // class HilbertSpaceHubbardHolstein

template<typename Word>
SizeType HilbertSpaceHubbardHolstein<Word>::numberBitphonons_ = 1;
} // namespace Dmrg
/*@}*/
#endif

