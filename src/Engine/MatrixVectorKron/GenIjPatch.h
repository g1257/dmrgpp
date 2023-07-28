/*
Copyright (c) 2012, UT-Battelle, LLC
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

/*! \file GenIjPatch.h
 *
 *
 */

#ifndef GEN_IJ_PATCH_HEADER_H
#define GEN_IJ_PATCH_HEADER_H

#include "Vector.h"
#include <cassert>

namespace Dmrg
{

template <typename LeftRightSuperType_>
class GenIjPatch
{

public:

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef typename LeftRightSuperType::BasisType BasisType;
	typedef typename BasisType::QnType QnType;
	typedef PsimagLite::Vector<SizeType>::Type VectorSizeType;

	enum LeftOrRightEnumType { LEFT = 0,
		RIGHT = 1 };

	GenIjPatch(const LeftRightSuperType& lrs, const QnType& target)
	    : lrs_(lrs)
	    , qn_(target)
	{
		for (SizeType i = 0; i < lrs.left().partition() - 1; i++) {
			for (SizeType j = 0; j < lrs.right().partition() - 1; ++j) {

				if (QnType(lrs.left().qnEx(i), lrs.right().qnEx(j)) != target)
					continue;

				patchesLeft_.push_back(i);
				patchesRight_.push_back(j);
			}
		}
	}

	const QnType& qn() const { return qn_; }

	const VectorSizeType& operator()(LeftOrRightEnumType leftOrRight) const
	{
		return (leftOrRight == LEFT) ? patchesLeft_ : patchesRight_;
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

private:

	const LeftRightSuperType& lrs_;
	const QnType& qn_;
	VectorSizeType patchesLeft_;
	VectorSizeType patchesRight_;

}; // class GenIjPatch
} // namespace PsimagLite

/*@}*/

#endif // GEN_IJ_PATCH_HEADER_H
