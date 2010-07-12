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

/*! \file ApplyOperatorLocal.h
 *
 *  documentation FIXME
 */
#ifndef APPLY_OPERATOR_LOCAL_H
#define APPLY_OPERATOR_LOCAL_H

#include "Utils.h"
#include "FermionSign.h"
#include "ProgramGlobals.h"

namespace Dmrg {
	
	template<typename BasisWithOperatorsType,typename VectorWithOffsetType_,typename TargetVectorType>
	class ApplyOperatorLocal {
			
			
			typedef typename BasisWithOperatorsType::RealType RealType;

		public:
			typedef typename BasisWithOperatorsType::BasisType BasisType;
			typedef VectorWithOffsetType_ VectorWithOffsetType;
			typedef typename BasisWithOperatorsType::OperatorType OperatorType;
			
			ApplyOperatorLocal(
					const BasisType& basisS,
       					const BasisType& basisE,
	    				const BasisType& basisSE)
			: basisS_(basisS),basisE_(basisE),basisSE_(basisSE)
			{
				//std::cerr<<"APPLY "<<basisSE_.size()<<" "<<basisSE.size()<<"\n";
			}
					//! FIXME: we need to make a fast version for when we're just
			//! figuring out where the (non-zero) partition is
			void operator()(
					VectorWithOffsetType& dest,
					const VectorWithOffsetType& src,
				     	const OperatorType& A,
	  				const FermionSign& fermionSign,
					size_t systemOrEnviron) const
			{
				//std::cerr<<"APPLY-OPERTOR() "<<basisSE_.size()<<"\n";
				if (systemOrEnviron == ProgramGlobals::EXPAND_SYSTEM) applyLocalOpSystem(dest,src,A,fermionSign);
				else applyLocalOpEnviron(dest,src,A);
			}
			
		private:
			void applyLocalOpSystem(
					VectorWithOffsetType& dest,
					const VectorWithOffsetType& src,
					const OperatorType& A,
	  				const FermionSign& fermionSign) const
			{
				TargetVectorType dest2(basisSE_.size());
				
				for (size_t i=0;i<dest2.size();i++) dest2[i] = 0;
				
				for (size_t ii=0;ii<src.sectors();ii++) {
					size_t i = src.sector(ii);
					applyLocalOpSystem(dest2,src,A,fermionSign,i);
				}
				dest.fromFull(dest2,basisSE_);
			}

			void applyLocalOpSystem(
					TargetVectorType& dest2,
					const VectorWithOffsetType& src,
					const OperatorType& A,
	  				const FermionSign& fermionSign,
					size_t i0) const
			{
				size_t offset = src.offset(i0);
				size_t final = offset + src.effectiveSize(i0);
				//size_t counter=0;
				size_t ns = basisS_.permutationVector().size();
				size_t nx = ns/A.data.rank();
				if (src.size()!=basisSE_.permutationVector().size()) throw std::runtime_error("applyLocalOpSystem SE\n");
				
				for (size_t i=offset;i<final;i++) {
					size_t x=0,y=0;
					utils::getCoordinates(x,y,basisSE_.permutation(i),ns);
					//if (y>=basisE_.permutationVector().size()) throw std::runtime_error("applyLocalOpSystem E\n");
					size_t x0=0,x1=0;
					if (x>=basisS_.permutationVector().size()) throw std::runtime_error("applyLocalOpSystem S\n");
					utils::getCoordinates(x0,x1,basisS_.permutation(x),nx);
					/*int nx0 = basisS_.electrons(x)-electrons[x1];
					if (nx0<0) throw std::runtime_error("TimeStepTargetting::applyLocalOpSystem(...)\n");
					*/
					RealType sign = fermionSign(x,A.fermionSign);
					for (int k=A.data.getRowPtr(x1);k<A.data.getRowPtr(x1+1);k++) {
						size_t x1prime = A.data.getCol(k);
						size_t xprime = basisS_.permutationInverse(x0+x1prime*nx);
						size_t j = basisSE_.permutationInverse(xprime+y*ns);
						dest2[j] += src[i]*A.data.getValue(k)*sign;
					}
				}
				
			}

			void applyLocalOpEnviron(
					VectorWithOffsetType& dest,
					const VectorWithOffsetType& src,
					const OperatorType& A) const
			{
				TargetVectorType dest2(basisSE_.size());
				
				for (size_t i=0;i<dest2.size();i++) dest2[i] = 0;
				
				for (size_t ii=0;ii<src.sectors();ii++) {
					size_t i = src.sector(ii);
					applyLocalOpEnviron(dest2,src,A,i);
				}
				dest.fromFull(dest2,basisSE_);
			}

			void applyLocalOpEnviron(
					TargetVectorType& dest2,
					const VectorWithOffsetType& src,
				     	const OperatorType& A,
					size_t i0) const
			{
				size_t offset = src.offset(i0);
				size_t final = offset + src.effectiveSize(i0);
				
				size_t ns = basisS_.size();
				size_t nx = A.data.rank();
				
				for (size_t i=offset;i<final;i++) {
					size_t x=0,y=0;
					utils::getCoordinates(x,y,basisSE_.permutation(i),ns);
					size_t y0=0,y1=0;
					utils::getCoordinates(y0,y1,basisE_.permutation(y),nx);
					RealType sign = basisS_.fermionicSign(x,A.fermionSign);
					for (int k=A.data.getRowPtr(y0);k<A.data.getRowPtr(y0+1);k++) {
						size_t y0prime = A.data.getCol(k);
						size_t yprime = basisE_.permutationInverse(y0prime+y1*nx);
						size_t j = basisSE_.permutationInverse(x+yprime*ns);
						dest2[j] += src[i]*A.data.getValue(k)*sign;
					}
				}
			}

			const BasisType& basisS_;
			const BasisType& basisE_;
			const BasisType& basisSE_;
	}; // class ApplyOperatorLocal
} // namespace Dmrg 

/*@}*/
#endif //APPLY_OPERATOR_LOCAL_H
