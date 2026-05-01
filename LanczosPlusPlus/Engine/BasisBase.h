
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2014 , UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef LANCZOS_BASIS_BASE_H
#define LANCZOS_BASIS_BASE_H
#include "LabeledOperator.h"
#include "LanczosGlobals.h"
#include "Vector.h"

namespace LanczosPlusPlus {

template <typename GeometryType> class BasisBase {

public:

	enum PrintEnum
	{
		PRINT_BINARY,
		PRINT_DECIMAL
	};

	typedef LanczosGlobals::PairIntType                 PairIntType;
	typedef LanczosGlobals::WordType                    WordType;
	typedef typename PsimagLite::Vector<WordType>::Type VectorWordType;
	typedef LabeledOperator                             LabeledOperatorType;

	virtual ~BasisBase() { }

	virtual PairIntType parts() const = 0;

	virtual SizeType dofs() const = 0;

	virtual SizeType size() const = 0;

	virtual SizeType hilbertOneSite(SizeType = 0) const = 0;

	virtual WordType operator()(SizeType i, SizeType spin) const = 0;

	virtual SizeType perfectIndex(const PsimagLite::Vector<WordType>::Type& kets) const = 0;

	virtual SizeType perfectIndex(WordType ket1, WordType ket2) const = 0;

	virtual SizeType perfectIndex(WordType newKet, SizeType ispace, SizeType spinOfNew) const
	    = 0;

	virtual PairIntType getBraIndex(WordType ket1,
	                                WordType ket2,
	                                const LabeledOperatorType&,
	                                SizeType site,
	                                SizeType spin,
	                                SizeType orb) const
	    = 0;

	virtual int doSign(WordType ket1,
	                   WordType ket2,
	                   SizeType i,
	                   SizeType orb1,
	                   SizeType j,
	                   SizeType orb2,
	                   SizeType spin) const
	    = 0;

	virtual int
	doSignGf(WordType a, WordType b, SizeType ind, SizeType spin, SizeType orb) const
	    = 0;

	virtual int
	doSignSpSm(WordType a, WordType b, SizeType ind, SizeType spin, SizeType orb) const
	{
		return 1;
	}

	virtual SizeType isThereAnElectronAt(WordType ket1,
	                                     WordType ket2,
	                                     SizeType site,
	                                     SizeType spin,
	                                     SizeType orb) const
	    = 0;

	virtual SizeType orbsPerSite(SizeType i) const = 0;

	virtual SizeType orbs() const = 0;

	virtual SizeType
	getN(WordType ket1, WordType ket2, SizeType site, SizeType spin, SizeType orb) const
	    = 0;

	virtual bool
	getBra(WordType&, WordType, WordType, const LabeledOperatorType&, SizeType, SizeType) const
	    = 0;

	virtual void print(std::ostream&, PrintEnum) const = 0;
}; // class BasisBase

} // namespace LanczosPlusPlus
#endif
