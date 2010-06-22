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

/*! \file FreeOperator.h
 *
 * 
 *
 */
#ifndef FREE_OPERATOR_H
#define FREE_OPERATOR_H

#include "Utils.h"

namespace Dmrg {
	// All interactions == 0
	template<typename MatrixType,typename HilbertVectorType>
	class FreeOperator {
			typedef unsigned int long long UnsignedIntegerType;
			//static size_t const SPIN_UP=0,SPIN_DOWN=1;
			typedef FlavoredState<UnsignedIntegerType> FlavoredStateType;
			typedef typename MatrixType::value_type FieldType;
			typedef typename HilbertVectorType::HilbertTermType HilbertTermType;
		public:
			FreeOperator(const MatrixType& U,const std::string& label,size_t site,size_t flavor,size_t dof) :
				U_(U),label_(label),site_(site),flavor_(flavor),dof_(dof)
			{
			}
			
			void apply(
				
				HilbertVectorType& dest,
				const HilbertVectorType& src) const
			{
				for (size_t i=0;i<src.terms();i++) {
					apply(dest,src.term(i));
				}
			}
			
			void apply(
					HilbertVectorType& dest,
					const HilbertTermType& src) const
			{
				dest.clear();
				size_t size = U_.n_row()/dof_;
				typename HilbertVectorType::HilbertTermType term = src;
				
				for (size_t lambda = 0;lambda < U_.n_col();lambda++) {
					
					applyInternal(term,src,lambda); // term.value contains sign
					FieldType factor = U_(site_+flavor_*size,lambda);
					term.value *= factor;
					dest.add(term);
				}
			}

		private:
			template<typename HilbertTermType>
			void applyInternal(
					HilbertTermType& dest,
					const HilbertTermType& src,
					size_t lambda) const
			{
				dest = src;
				int sign = dest.state.apply(label_,flavor_,lambda);
				dest.value *= sign; 
			}
			
			const MatrixType& U_;
			std::string label_;
			size_t site_;
			size_t flavor_;
			size_t dof_;
	}; // FreeOperator
} // namespace Dmrg 

/*@}*/
#endif
