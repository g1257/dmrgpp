// BEGIN LICENSE BLOCK
/*
Copyright (c) 2008 , UT-Battelle, LLC
All rights reserved

[DMRG++, Version 1.0.0]
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
#include "Profiling.h"

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

		enum {LEFT_BRACKET=ObserverHelperType::LEFT_BRACKET,
			RIGHT_BRACKET=ObserverHelperType::RIGHT_BRACKET};
	public:
		OnePointCorrelations(ObserverHelperType& helper,
				bool verbose=false)
		: helper_(helper),
		  verbose_(verbose)
		{}

		template<typename ApplyOperatorType>
		FieldType operator()(size_t site,
			const typename ApplyOperatorType::OperatorType& A,
			bool corner = false)
		{
			size_t pnter=site;
			helper_.setPointer(pnter);
			try {
				const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET);
				const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET);

				return onePointInternal<ApplyOperatorType>(site,A,src1,src2,corner);
			} catch (std::exception& e) {
				std::cerr<<"CAUGHT: "<<e.what();
				std::cerr<<"WARNING: Observer::onePoint(...): Nothing here yet\n";
				return 0;
			}
		}

		template<typename ApplyOperatorType>
		FieldType hookForZero(size_t site,
				      const typename ApplyOperatorType::OperatorType& A,
				      bool corner = false)
		{
			size_t pnter=site;
			helper_.setPointer(pnter);
			try {
				const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(LEFT_BRACKET);
				const VectorWithOffsetType& src2 =  helper_.getVectorFromBracketId(RIGHT_BRACKET);

				return onePointInternalHookForZero<ApplyOperatorType>(site,A,src1,src2,corner);
			} catch (std::exception& e) {
				std::cerr<<"CAUGHT: "<<e.what();
				std::cerr<<"WARNING: Observer::onePoint(...): Nothing here yet\n";
				return 0;
			}
		}

	private:

		template<typename ApplyOperatorType>
		FieldType onePointInternal(size_t site,
					   const typename ApplyOperatorType::OperatorType& A,
					   const VectorWithOffsetType& src1,
					   const VectorWithOffsetType& src2,
					   bool corner = false)
		{
			
			ApplyOperatorType applyOpLocal1(helper_.leftRightSuper());
			VectorWithOffsetType dest;
//			assert(helper_.fermionicSignLeft().size()==helper_.leftRightSuper().left().size());
			applyOpLocal1(dest,src1,A,helper_.fermionicSignLeft(),
					helper_.direction(),corner);
				
			FieldType sum = static_cast<FieldType>(0.0);
			const VectorWithOffsetType& v1 = dest;
			const VectorWithOffsetType& v2 = src2;
			for (size_t ii=0;ii<v1.sectors();ii++) {
				size_t i = v1.sector(ii);
				for (size_t jj=0;jj<v1.sectors();jj++) {
					size_t j = v2.sector(jj);
					if (i!=j) continue;
					size_t offset = v1.offset(i);
					for (size_t k=0;k<v1.effectiveSize(i);k++) 
						sum+= v1[k+offset] * std::conj(v2[k+offset]);
				}
			}
			return sum;
		}

		template<typename ApplyOperatorType>
		FieldType onePointInternalHookForZero(size_t site,
						      const typename ApplyOperatorType::OperatorType& A,
						      const VectorWithOffsetType& src1,
						      const VectorWithOffsetType& src2,
						      bool corner = false)
		{

			ApplyOperatorType applyOpLocal1(helper_.leftRightSuper());
			VectorWithOffsetType dest;
//			assert(helper_.fermionicSignLeft().size()==helper_.leftRightSuper().left().size());
			applyOpLocal1.hookForZero(dest,src1,A,helper_.fermionicSignLeft(),helper_.direction(),corner);

			FieldType sum = static_cast<FieldType>(0.0);
			const VectorWithOffsetType& v1 = dest;
			const VectorWithOffsetType& v2 = src2;
			for (size_t ii=0;ii<v1.sectors();ii++) {
				size_t i = v1.sector(ii);
				for (size_t jj=0;jj<v1.sectors();jj++) {
					size_t j = v2.sector(jj);
					if (i!=j) continue;
					size_t offset = v1.offset(i);
					for (size_t k=0;k<v1.effectiveSize(i);k++)
						sum+= v1[k+offset] * std::conj(v2[k+offset]);
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
