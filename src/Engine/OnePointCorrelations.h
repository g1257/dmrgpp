/*
Copyright (c) 2008-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file OnePointCorrelations.h
 *
 *
 *  A class to perform post-processing calculation of one-point correlations
 *  <state1 | A_i | state2>
 * Note that site below is not necessarily the site, but there's a mapping
 * mapping(site) = i
 */

#ifndef ONE_POINT_H
#define ONE_POINT_H
#include "CrsMatrix.h"
#include "GetBraOrKet.h"
#include "VectorWithOffset.h" // for operator*
#include "VectorWithOffsets.h" // for operator*

namespace Dmrg
{

template <typename ObserverHelperType, typename ModelType>
class OnePointCorrelations
{

	typedef typename ObserverHelperType::MatrixType MatrixType;
	typedef typename ObserverHelperType::VectorType VectorType;
	typedef typename ObserverHelperType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename VectorType::value_type FieldType;
	typedef typename BasisWithOperatorsType::RealType RealType;

public:

	OnePointCorrelations(const ObserverHelperType& helper,
	    const ModelType& model)
	    : helper_(helper)
	    , model_(model)
	{
	}

	template <typename ApplyOperatorType>
	FieldType operator()(SizeType ptr,
	    const typename ApplyOperatorType::OperatorType& A,
	    SizeType site,
	    typename ApplyOperatorType::BorderEnum corner,
	    const PsimagLite::GetBraOrKet& bra,
	    const PsimagLite::GetBraOrKet& ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);

			return onePointInternal<ApplyOperatorType>(A, site, src1, src2, corner, ptr);
		} catch (std::exception& e) {
			std::cerr << "CAUGHT: " << e.what();
			std::cerr << "WARNING: Observer::onePoint(...): Nothing here yet\n";
			return 0;
		}
	}

	template <typename ApplyOperatorType>
	FieldType hookForZero(SizeType ptr,
	    const typename ApplyOperatorType::OperatorType& A,
	    SizeType splitSize,
	    const PsimagLite::GetBraOrKet& bra,
	    const PsimagLite::GetBraOrKet& ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);

			return onePointInternalHookForZero<ApplyOperatorType>(A,
			    splitSize,
			    src1,
			    src2,
			    ptr);
		} catch (std::exception& e) {
			std::cerr << "CAUGHT: " << e.what();
			std::cerr << "WARNING: Observer::onePoint(...): Nothing here yet\n";
			return 0;
		}
	}

private:

	template <typename ApplyOperatorType>
	FieldType onePointInternal(const typename ApplyOperatorType::OperatorType& A,
	    SizeType site,
	    const VectorWithOffsetType& src1,
	    const VectorWithOffsetType& src2,
	    typename ApplyOperatorType::BorderEnum corner,
	    SizeType ptr) const
	{
		if (src1.sectors() == 0 || src2.sectors() == 0)
			return 0.0;

		SizeType splitSize = model_.hilbertSize(site);
		ApplyOperatorType applyOpLocal1(helper_.leftRightSuper(ptr),
		    true);
		VectorWithOffsetType dest;
		applyOpLocal1(dest,
		    src1,
		    A,
		    helper_.fermionicSignLeft(ptr),
		    splitSize,
		    helper_.direction(ptr),
		    corner);

		FieldType sum = static_cast<FieldType>(0.0);
		const VectorWithOffsetType& v1 = dest;
		const VectorWithOffsetType& v2 = src2;
		for (SizeType ii = 0; ii < v1.sectors(); ii++) {
			SizeType i = v1.sector(ii);
			for (SizeType jj = 0; jj < v1.sectors(); jj++) {
				SizeType j = v2.sector(jj);
				if (i != j)
					continue;
				SizeType offset = v1.offset(i);
				for (SizeType k = 0; k < v1.effectiveSize(i); k++)
					sum += v1.slowAccess(k + offset) * PsimagLite::conj(v2.slowAccess(k + offset));
			}
		}

		return sum;
	}

	template <typename ApplyOperatorType>
	FieldType onePointInternalHookForZero(const typename ApplyOperatorType::OperatorType& A,
	    SizeType splitSize,
	    const VectorWithOffsetType& src1,
	    const VectorWithOffsetType& src2,
	    SizeType ptr) const
	{

		ApplyOperatorType applyOpLocal1(helper_.leftRightSuper(ptr),
		    true);
		VectorWithOffsetType dest;
		applyOpLocal1.hookForZero(dest,
		    src1,
		    A,
		    splitSize,
		    // helper_.fermionicSignLeft(ptr),
		    helper_.direction(ptr));

		FieldType sum = static_cast<FieldType>(0.0);
		const VectorWithOffsetType& v1 = dest;
		const VectorWithOffsetType& v2 = src2;
		for (SizeType ii = 0; ii < v1.sectors(); ii++) {
			SizeType i = v1.sector(ii);
			for (SizeType jj = 0; jj < v1.sectors(); jj++) {
				SizeType j = v2.sector(jj);
				if (i != j)
					continue;
				SizeType offset = v1.offset(i);
				for (SizeType k = 0; k < v1.effectiveSize(i); k++)
					sum += v1.slowAccess(k + offset) * PsimagLite::conj(v2.slowAccess(k + offset));
			}
		}

		return sum;
	}

	const ObserverHelperType& helper_;
	const ModelType& model_;
}; // class OnePointCorrelations
} // namespace Dmrg

/*@}*/
#endif // ONE_POINT_H
