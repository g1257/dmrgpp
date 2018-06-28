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
#include "VectorWithOffsets.h" // for operator*
#include "VectorWithOffset.h" // for operator*

namespace Dmrg {

template<typename ObserverHelperType>
class OnePointCorrelations {
	typedef typename ObserverHelperType::MatrixType MatrixType;
	typedef typename ObserverHelperType::VectorType VectorType ;
	typedef typename ObserverHelperType::VectorWithOffsetType
	VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType
	BasisWithOperatorsType;
	typedef typename VectorType::value_type FieldType;
	typedef typename BasisWithOperatorsType::RealType RealType;

	enum {LEFT_BRAKET=ObserverHelperType::LEFT_BRAKET,
		  RIGHT_BRAKET=ObserverHelperType::RIGHT_BRAKET};
public:
	OnePointCorrelations(ObserverHelperType& helper,
	                     bool verbose=false)
	    : helper_(helper),
	      verbose_(verbose)
	{}

	template<typename ApplyOperatorType>
	FieldType operator()(SizeType site,
	                     const typename ApplyOperatorType::OperatorType& A,
	                     typename ApplyOperatorType::BorderEnum corner)
	{
		SizeType threadId =0;
		SizeType pnter=site;
		helper_.setPointer(threadId,pnter);
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRAKET,
			                                                                  threadId);
			const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRAKET,
			                                                                   threadId);

			return onePointInternal<ApplyOperatorType>(site,A,src1,src2,corner,threadId);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: Observer::onePoint(...): Nothing here yet\n";
			return 0;
		}
	}

	template<typename ApplyOperatorType>
	FieldType hookForZero(SizeType site,
	                      const typename ApplyOperatorType::OperatorType& A,
	                      bool corner = false)
	{
		SizeType pnter=site;
		SizeType threadId = 0;
		helper_.setPointer(threadId,pnter);
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRAKET,
			                                                                  threadId);
			const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRAKET,
			                                                                   threadId);

			return onePointInternalHookForZero<ApplyOperatorType>(site,
			                                                      A,
			                                                      src1,
			                                                      src2,
			                                                      corner,
			                                                      threadId);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: Observer::onePoint(...): Nothing here yet\n";
			return 0;
		}
	}

private:

	template<typename ApplyOperatorType>
	FieldType onePointInternal(SizeType,
	                           const typename ApplyOperatorType::OperatorType& A,
	                           const VectorWithOffsetType& src1,
	                           const VectorWithOffsetType& src2,
	                           typename ApplyOperatorType::BorderEnum corner,
	                           SizeType threadId)
	{
		if (src1.sectors() == 0 || src2.sectors() == 0) return 0.0;
		ApplyOperatorType applyOpLocal1(helper_.leftRightSuper(threadId),
		                                helper_.withLegacyBugs());
		VectorWithOffsetType dest;
		applyOpLocal1(dest,
		              src1,
		              A,
		              helper_.fermionicSignLeft(threadId),
		              helper_.direction(threadId),corner);

		FieldType sum = static_cast<FieldType>(0.0);
		const VectorWithOffsetType& v1 = dest;
		const VectorWithOffsetType& v2 = src2;
		for (SizeType ii=0;ii<v1.sectors();ii++) {
			SizeType i = v1.sector(ii);
			for (SizeType jj=0;jj<v1.sectors();jj++) {
				SizeType j = v2.sector(jj);
				if (i!=j) continue;
				SizeType offset = v1.offset(i);
				for (SizeType k=0;k<v1.effectiveSize(i);k++)
					sum+= v1.slowAccess(k+offset)*PsimagLite::conj(v2.slowAccess(k+offset));
			}
		}
		return sum;
	}

	template<typename ApplyOperatorType>
	FieldType onePointInternalHookForZero(SizeType,
	                                      const typename ApplyOperatorType::OperatorType& A,
	                                      const VectorWithOffsetType& src1,
	                                      const VectorWithOffsetType& src2,
	                                      bool, //= false
	                                      SizeType threadId)
	{

		ApplyOperatorType applyOpLocal1(helper_.leftRightSuper(threadId),
		                                helper_.withLegacyBugs());
		VectorWithOffsetType dest;
		applyOpLocal1.hookForZero(dest,
		                          src1,
		                          A,
		                          helper_.fermionicSignLeft(threadId),
		                          helper_.direction(threadId));

		FieldType sum = static_cast<FieldType>(0.0);
		const VectorWithOffsetType& v1 = dest;
		const VectorWithOffsetType& v2 = src2;
		for (SizeType ii=0;ii<v1.sectors();ii++) {
			SizeType i = v1.sector(ii);
			for (SizeType jj=0;jj<v1.sectors();jj++) {
				SizeType j = v2.sector(jj);
				if (i!=j) continue;
				SizeType offset = v1.offset(i);
				for (SizeType k=0;k<v1.effectiveSize(i);k++)
					sum+= v1.slowAccess(k+offset)*PsimagLite::conj(v2.slowAccess(k+offset));
			}
		}
		return sum;
	}

	ObserverHelperType& helper_;
	bool verbose_;

};  //class OnePointCorrelations
} // namespace Dmrg

/*@}*/
#endif // ONE_POINT_H
