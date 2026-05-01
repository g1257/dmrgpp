/*
Copyright (c) 2009-2012, UT-Battelle, LLC
All rights reserved

[LanczosPlusPlus, Version 1.0.0]
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

/** \ingroup LanczosPlusPlus */
/*@{*/

/*! \file LanczosGlobals.h
 *
 *
 *
 */
#ifndef LANCZOS_PROGRAM_LIMITS_H
#define LANCZOS_PROGRAM_LIMITS_H
#include "BitManip.h"
#include "Combinatorial.hh"
#include "CrsMatrix.h"
#include "TypeToString.h"
#include "Vector.h"
#include <bitset>
#include <climits>

namespace LanczosPlusPlus {

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1, T2>& p)
{
	os << p.first << " " << p.second;
	return os;
}

struct LanczosGlobals {

	typedef std::pair<int, int> PairIntType;
	typedef unsigned int long   WordType;

	static int FERMION_SIGN;

	enum
	{
		FERMION,
		BOSON
	};

	enum
	{
		SPIN_UP,
		SPIN_DOWN
	};

	enum ConnectionEnum
	{
		NONE
	};

	static int doSign(WordType a, SizeType i)
	{
		WordType mask = (1 << i) - 1;
		// Parity of single occupied between i and nsite-1
		return (PsimagLite::BitManip::count(a & mask) & 1) ? FERMION_SIGN : 1;
	}

	template <typename T> static void binRep(std::ostream& os, SizeType n, const T& a)
	{
		const char*        beg = reinterpret_cast<const char*>(&a);
		SizeType           len = sizeof(a);
		const char*        end = beg + len;
		PsimagLite::String buffer("");

		while (beg != end) {
			PsimagLite::String str(std::bitset<CHAR_BIT>(*(end - 1)).to_string());
			buffer += str;
			end--;
		}

		SizeType lmn = buffer.length();
		if (lmn >= n)
			lmn -= n;
		os << buffer.substr(lmn, n);
	}

	template <typename VectorWordType>
	static void printBasisBinary(std::ostream& os, SizeType n, VectorWordType& data)
	{
		for (SizeType i = 0; i < data.size(); i++) {
			binRep(os, n, data[i]);
			os << "\n";
		}

		os << "--------------\n";
	}

	static void printBasisBinary(std::ostream& os, SizeType n, SizeType h)
	{
		for (SizeType i = 0; i < h; ++i) {
			binRep(os, n, i);
			os << "\n";
		}

		os << "--------------\n";
	}

	template <typename WordType>
	static void
	printBasisDecimal(std::ostream& os, SizeType n, const std::vector<WordType>& data)
	{
		for (SizeType i = 0; i < data.size(); i++) {
			os << data[i] << " ";
			if (i > 0 && i % n == 0)
				std::cout << "\n";
		}

		os << "\n--------------\n";
	}

	template <typename WordType>
	static void printBasisDecimal(std::ostream&                                     os,
	                              SizeType                                          n,
	                              const std::vector<std::pair<WordType, WordType>>& data)
	{
		for (SizeType i = 0; i < data.size(); i++) {
			os << data[i].first << " " << data[i].second;
			if (i > 0 && i % n == 0)
				std::cout << "\n";
		}

		os << "\n--------------\n";
	}

	static void printBasisDecimal(std::ostream& os, SizeType n, SizeType h)
	{
		for (SizeType i = 0; i < h; ++i) {
			os << i << " ";
			if (i > 0 && i % n == 0)
				std::cout << "\n";
		}

		os << "\n--------------\n";
	}

	template <typename SomeVectorType>
	static typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
	                                     void>::Type
	transform(SomeVectorType&                                                   gs,
	          SizeType                                                          offset,
	          SomeVectorType&                                                   gstmp,
	          const PsimagLite::CrsMatrix<typename SomeVectorType::value_type>& tr)
	{
		for (SizeType i = 0; i < gs.size(); i++) {
			assert(i + offset < gstmp.size());
			gstmp[i + offset] = gs[i];
		}

		PsimagLite::CrsMatrix<typename SomeVectorType::value_type> rT;
		transposeConjugate(rT, tr);
		gs.clear();
		gs.resize(tr.rows());
		multiply(gs, rT, gstmp);
	}

	static void doBitmask(SizeType total)
	{
		if (total == bitmask_.size())
			return;

		bitmask_.resize(total);
		bitmask_[0] = 1ul;
		for (SizeType i = 1; i < bitmask_.size(); i++)
			bitmask_[i] = bitmask_[i - 1] << 1;
	}

	static const WordType& bitmask(SizeType i)
	{
		assert(i < bitmask_.size());
		return bitmask_[i];
	}

	static void doCombinatorial(SizeType total)
	{
		if (total == comb_.size()) {
			return;
		}

		comb_.resize(total);
	}

	static const SizeType& combinatorial(SizeType i, SizeType j) { return comb_(i, j); }

private:

	static PsimagLite::Vector<WordType>::Type bitmask_;
	static Combinatorial                      comb_;
}; // LanczosGlobals

} // namespace LanczosPlusPlus
/*@}*/
#endif
