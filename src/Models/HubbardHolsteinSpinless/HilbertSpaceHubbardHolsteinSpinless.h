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
 *  This class represents the Hilbert space for the Hubbard Holstein Model Spinless
 *  States are represented with binary numbers. N bits per site
 * (1 for the fermions, N-1 for the phonons)
 *  Bits meaning:
 *  0....00  empty state
 *  0..1.00  a 1 at location x>1 means state with N-1 maxphonons and zero fermions
 *  1111111  all ones means 1 electrons
 *
 *  Note: this is a static class
 *
 */
#ifndef HILBERTSPACE_HUBBARD_HOLSTEIN_SPINLESS_HEADER_H
#define HILBERTSPACE_HUBBARD_HOLSTEIN_SPINLESS_HEADER_H

#include "Utils.h"

namespace Dmrg {

//! A class to operate on n-ary numbers (base n)
template <typename Word> class HilbertSpaceHubbardHolsteinSpinless {

public:

	using HilbertState = Word;

	// Get electronic state from combined electron and phonon ket a
	static Word getF(Word a) { return (a & 1); }

	// Get phononic state from combined electron and phonon ket a
	static Word getP(Word a) { return (a >> 1); }

	// Create electron in binary number "a"
	static void createF(Word& a)
	{
		Word mask = 1;
		assert((a & mask) == 0);
		a |= mask;
	}

	// Create phonon
	static void createP(Word& a)
	{
		const SizeType nphonons = getP(a);
		const Word     stateP   = 1 + nphonons;
		const Word     maskP    = (stateP << 1);
		const Word     maskF    = getF(a);
		assert(maskF < 2);
		a = (maskP | maskF);
	}

	static SizeType electronsWithGivenSpin(Word a)
	{
		const SizeType f = getF(a);
		return (f & 1);
	}
}; // class HilbertSpaceHubbardHolstein

} // namespace Dmrg
/*@}*/
#endif
