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

#include <vector>
#include <cassert>
#include "GenGroup.h"
//#include "QvalStruct.h"

namespace Dmrg {

template<typename LeftRightSuperType>
class GenIjPatch {

	typedef typename LeftRightSuperType::BasisType BasisType;

public:

	typedef GenGroup<BasisType> GenGroupType;
//	typedef QvalStruct<LeftRightSuperType> QvalStructType;

	enum LeftOrRightEnumType {LEFT=0,RIGHT=1};

	GenIjPatch(const LeftRightSuperType& lrs,int target)
	{
		GenGroupType groupLeft(lrs.left());
		GenGroupType groupRight(lrs.right());

//		std::cerr<<"groupLeft.size="<<groupLeft.size()<<"\n";
//		std::cerr<<"groupRight.size="<<groupRight.size()<<"\n";

		for (size_t i=0;i<groupLeft.size()-1;i++) {
			size_t istart = groupLeft(i);
			assert(istart<lrs.left().size());
//			size_t iend = groupLeft(i+1)-1;
			for (size_t j=0;j<groupRight.size()-1;j++) {
				size_t jstart = groupRight(j);
				assert(jstart<lrs.right().size());
//				size_t jend = groupRight(j+1);

				if (lrs.left().qn(istart) + lrs.right().qn(jstart)!=target) continue;

//				npatches++;
				patchesLeft_.push_back(i);
				patchesRight_.push_back(j);
			}
		}
	}

	size_t operator()(LeftOrRightEnumType leftOrRight,size_t i) const
	{
		assert(i<patchesLeft_.size() && i<patchesRight_.size());
		return (leftOrRight==LEFT) ? patchesLeft_[i] : patchesRight_[i];
	}

	size_t size() const
	{
		assert(patchesLeft_.size()==patchesRight_.size());
		return patchesLeft_.size();
	}

private:

	std::vector<size_t> patchesLeft_,patchesRight_;

}; //class GenIjPatch
} // namespace PsimagLite

/*@}*/

#endif // GEN_IJ_PATCH_HEADER_H
