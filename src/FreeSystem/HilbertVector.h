// BEGIN LICENSE BLOCK
/*
Copyright © 2009 , UT-Battelle, LLC
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

/*! \file HilbertVector.h
 *
 * Raw computations for a free Hubbard model
 *
 */
#ifndef HILBERT_VECTOR_H
#define HILBERT_VECTOR_H

#include "Utils.h"
#include "FlavoredState.h"

namespace Dmrg {
	// All interactions == 0
	
	template<typename FieldType,typename FlavoredStateType>
	struct HilbertTerm {
		const FlavoredStateType* state;
		const FieldType* value;	
	};
	
	template<typename FieldType>
	class HilbertVector {
			typedef unsigned int long long UnsignedIntegerType;
			//static size_t const SPIN_UP=0,SPIN_DOWN=1;
			typedef FlavoredState<UnsignedIntegerType> FlavoredStateType;
			
			
		public:
			typedef HilbertTerm<FieldType,FlavoredStateType> HilbertTermType;
			
			HilbertVector(size_t size,size_t dof) :
				size_(size),dof_(dof)
			{
			}
			
			void add(const HilbertVector& another)
			{
				for (size_t i=0;i<another.data_.size();i++) {
					internalAdd(another.data_[i],another.values_[i]);
				}
			}
			
			HilbertTermType term(size_t i)
			{
				HilbertTermType t;
				t.state = &(data_[i]);
				t.value = &(values_[i]);
				return t;
			}

		private:
			
			void internalAdd(const FlavoredStateType& state,const FieldType& value)
			{
				size_t x = utils::isInVector(data_,state);
				if (x<0) {
					data_.push_back(state);
					values_.push_back(value);
				} else {
					values_[x] += value;
				}
			}
			
			size_t size_;
			size_t dof_;
			std::vector<FlavoredStateType> data_;
			std::vector<FieldType> values_;
			
	}; // HilbertVector
} // namespace Dmrg 

/*@}*/
#endif
