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

/*! \file HilbertSpaceFeAs.h
 *
 *  This class represents the Hilbert space for the FeAs Model
 *  States are represented with binary numbers. N bits per site
 *  Bits meaning:
 *  0.....0  empty state
 *  0..1..0  a 1 at location x means state "x"
 *  0..010..010..0, a 1 at location x and a 1 at location y means
 *  2 electrons on the site with states x and y respectively
 *  ...
 *  1111...111111  all ones means N electrons each with a different state
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACEFEAS_HEADER_H
#define HILBERTSPACEFEAS_HEADER_H

namespace Dmrg
{

//! A class to operate on n-ary numbers (base n)
template <typename Word>
class HilbertSpaceFeAs
{

	static SizeType orbitals_;

public:

	typedef Word HilbertState;

	enum { SPIN_UP = 0,
		SPIN_DOWN = 1 };

	static void setOrbitals(SizeType orbitals)
	{
		orbitals_ = orbitals;
	}

	// Get electronic state on site "j" in binary number "a"
	static Word get(Word const& a, SizeType j)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType k = dofs * j;
		SizeType ones = (1 << (dofs)) - 1;
		Word mask = (ones << k);

		mask &= a;
		mask >>= k;
		return mask;
	}

	// Create electron with internal dof "sigma" on site "j" in binary number "a"
	static void create(Word& a, SizeType j, SizeType sigma)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType k = dofs * j;
		Word mask = (1 << (k + sigma));
		a |= mask;
	}

	// Destroy electron with internal dof "sigma" on site "j" in binary number "a"
	static void destroy(Word& a, SizeType j, SizeType sigma)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType k = dofs * j;
		Word mask = (1 << (k + sigma));
		a &= (~mask);
	}

	// Is there an electron with internal dof "sigma" on site "i" in binary number "ket"?
	static bool isNonZero(Word const& ket, SizeType i, SizeType sigma)
	{

		Word tmp = get(ket, i);
		if (tmp & (1 << sigma))
			return true;

		return false;
	}

	//! returns the number of electrons of internal dof "value" in binary number "data"
	static int getNofDigits(Word const& data, SizeType value)
	{
		SizeType dofs = 2 * orbitals_;
		int ret = 0;
		Word data2 = data;
		SizeType i = 0;
		do {
			if ((data & (1 << (dofs * i + value))))
				ret++;
			i++;
		} while (data2 >>= dofs);

		return ret;
	}

	//! Number of electrons with spin spin (sums over bands and sites)
	static int electronsWithGivenSpin(Word const& data, SizeType spin)
	{
		SizeType sum = 0;
		Word data2 = data;
		SizeType digit = 0;
		while (data2 > 0) {
			SizeType sigma = digit % (2 * orbitals_);
			SizeType spin2 = static_cast<SizeType>(sigma / orbitals_);
			Word thisbit = (data2 & 1);
			data2 >>= 1;
			digit++;
			if (spin == spin2)
				sum += thisbit;
		}

		return sum;
	}

	//! Number of electrons in binary number "data" (sum over all bands)
	static SizeType electrons(const Word& data)
	{
		SizeType sum = 0;
		SizeType dofs = 2 * orbitals_;
		for (SizeType sector = 0; sector < dofs; sector++)
			sum += calcNofElectrons(data, sector);

		return sum;
	}
	//! Number of electrons with dof sector between i and j excluding
	//! i and j in binary number "ket"
	//!  intended for when i<j
	static int calcNofElectrons(Word const& ket, SizeType i, SizeType j, SizeType sector)
	{
		SizeType dofs = 2 * orbitals_;
		SizeType ii = i + 1;
		if (ii >= j)
			return 0;
		Word m = 0;
		for (SizeType k = dofs * ii; k < dofs * j; k++)
			m |= (1 << k);
		m = m & ket;
		return getNofDigits(m, sector);
	}

	//! Number of electrons with dof sector in binary number "ket"
	static int calcNofElectrons(Word const& ket, SizeType sector)
	{
		SizeType dofs = 2 * orbitals_;
		Word ket2 = ket;
		SizeType digit = 0;
		int sum = 0;
		while (ket2 > 0) {
			SizeType sector2 = digit % dofs;
			SizeType thisbit = (ket2 & 1);
			if (sector == sector2)
				sum += thisbit;
			digit++;
			ket2 >>= 1;
		}

		return sum;
	}

}; // class HilbertSpaceFeAs

template <typename Word>
SizeType HilbertSpaceFeAs<Word>::orbitals_ = 2;
} // namespace Dmrg

/*@}*/
#endif
