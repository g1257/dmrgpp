// BEGIN LICENSE BLOCK
/*
Copyright (c) 2009, UT-Battelle, LLC
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

#include "PackIndices.h" // in PsimagLite
#include "FermionSign.h"
#include "ProgramGlobals.h"

namespace Dmrg {
	
	template<typename LeftRightSuperType,typename VectorWithOffsetType_,typename TargetVectorType>
	class ApplyOperatorLocal {

		typedef typename LeftRightSuperType::BasisWithOperatorsType
				BasisWithOperatorsType;
		typedef typename BasisWithOperatorsType::RealType RealType;
		typedef PsimagLite::PackIndices PackIndicesType;

	public:
		enum {MIDDLE,LEFT_CORNER,RIGHT_CORNER};

		typedef typename BasisWithOperatorsType::BasisType BasisType;
		typedef VectorWithOffsetType_ VectorWithOffsetType;
		typedef typename BasisWithOperatorsType::OperatorType OperatorType;

		ApplyOperatorLocal(const LeftRightSuperType& lrs)
		: lrs_(lrs)
		{}

		//! FIXME: we need to make a fast version for when we're just
		//! figuring out where the (non-zero) partition is
		void operator()(
				VectorWithOffsetType& dest,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				const FermionSign& fermionSign,
				size_t systemOrEnviron,bool corner = false) const
		{
			if (!corner) {
				if (systemOrEnviron == ProgramGlobals::EXPAND_SYSTEM) applyLocalOpSystem(dest,src,A,fermionSign);
				else applyLocalOpEnviron(dest,src,A);
				return;
			}
			applyLocalOpCorner(dest,src,A,fermionSign);
		}

	private:
		void applyLocalOpSystem(
				VectorWithOffsetType& dest,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				const FermionSign& fermionSign,
				size_t whichPartOfTheLattice = MIDDLE) const
		{
			TargetVectorType dest2(lrs_.super().size());

			for (size_t i=0;i<dest2.size();i++) dest2[i] = 0;

			for (size_t ii=0;ii<src.sectors();ii++) {
				size_t i = src.sector(ii);
				switch (whichPartOfTheLattice) {
				case MIDDLE:
					applyLocalOpSystem(dest2,src,A,fermionSign,i);
					break;
				case LEFT_CORNER:
					throw std::runtime_error("applyLocalOpSystem: internal error\n");
					break;
				case RIGHT_CORNER:
					applyLocalOpRightCorner(dest2,src,A,i);
					break;
				}
			}
			dest.fromFull(dest2,lrs_.super());
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
			size_t ns = lrs_.left().permutationVector().size();
			size_t nx = ns/A.data.rank();
			if (src.size()!=lrs_.super().permutationVector().size())
				throw std::runtime_error("applyLocalOpSystem SE\n");

			PackIndicesType pack1(ns);
			PackIndicesType pack2(nx);
			for (size_t i=offset;i<final;i++) {
				size_t x=0,y=0;
				pack1.unpack(x,y,lrs_.super().permutation(i));
				//if (y>=basisE_.permutationVector().size()) throw std::runtime_error("applyLocalOpSystem E\n");
				size_t x0=0,x1=0;
				if (x>=lrs_.left().permutationVector().size())
					throw std::runtime_error("applyLocalOpSystem S\n");
				pack2.unpack(x0,x1,lrs_.left().permutation(x));
				/*int nx0 = basisS_.electrons(x)-electrons[x1];
					if (nx0<0) throw std::runtime_error("TimeStepTargetting::applyLocalOpSystem(...)\n");
				 */
				RealType sign = fermionSign(x,A.fermionSign);
				for (int k=A.data.getRowPtr(x1);k<A.data.getRowPtr(x1+1);k++) {
					size_t x1prime = A.data.getCol(k);
					size_t xprime = lrs_.left().permutationInverse(x0+x1prime*nx);
					size_t j = lrs_.super().permutationInverse(xprime+y*ns);
					dest2[j] += src[i]*A.data.getValue(k)*sign;
				}
			}

		}

		void applyLocalOpEnviron(
				VectorWithOffsetType& dest,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				size_t whichPartOfTheLattice = MIDDLE) const
		{
			TargetVectorType dest2(lrs_.super().size());

			for (size_t i=0;i<dest2.size();i++) dest2[i] = 0;

			for (size_t ii=0;ii<src.sectors();ii++) {
				size_t i = src.sector(ii);
				switch (whichPartOfTheLattice) {
				case MIDDLE:
					applyLocalOpEnviron(dest2,src,A,i);
					break;
				case LEFT_CORNER:
					applyLocalOpLeftCorner(dest2,src,A,i);
					break;
				case RIGHT_CORNER:
					throw std::runtime_error("applyLocalOpEnviron: internal error\n");
					break;
				}
			}
			dest.fromFull(dest2,lrs_.super());
		}

		void applyLocalOpEnviron(
				TargetVectorType& dest2,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				size_t i0) const
		{
			size_t offset = src.offset(i0);
			size_t final = offset + src.effectiveSize(i0);

			size_t ns = lrs_.left().size();
			size_t nx = A.data.rank();
			PackIndicesType pack1(ns);
			PackIndicesType pack2(nx);

			for (size_t i=offset;i<final;i++) {
				size_t x=0,y=0;
				pack1.unpack(x,y,lrs_.super().permutation(i));
				size_t y0=0,y1=0;
				pack2.unpack(y0,y1,lrs_.right().permutation(y));
				RealType sign = lrs_.left().fermionicSign(x,A.fermionSign);
				for (int k=A.data.getRowPtr(y0);k<A.data.getRowPtr(y0+1);k++) {
					size_t y0prime = A.data.getCol(k);
					size_t yprime = lrs_.right().permutationInverse(y0prime+y1*nx);
					size_t j = lrs_.super().permutationInverse(x+yprime*ns);
					dest2[j] += src[i]*A.data.getValue(k)*sign;
				}
			}
		}

		void applyLocalOpLeftCorner(
				TargetVectorType& dest2,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				size_t i0) const
		{
			size_t offset = src.offset(i0);
			size_t final = offset + src.effectiveSize(i0);

			size_t ns = lrs_.left().size();
			PackIndicesType pack(ns);

			for (size_t i=offset;i<final;i++) {
				size_t x=0,y=0;
				pack.unpack(x,y,lrs_.super().permutation(i));

				for (int k=A.data.getRowPtr(x);k<A.data.getRowPtr(x+1);k++) {
					size_t xprime = A.data.getCol(k);
					size_t j = lrs_.super().permutationInverse(xprime+y*ns);
					dest2[j] += src[i]*A.data.getValue(k);
				}
			}
		}

		void applyLocalOpRightCorner(
				TargetVectorType& dest2,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				size_t i0) const
		{
			size_t offset = src.offset(i0);
			size_t final = offset + src.effectiveSize(i0);
			//size_t counter=0;
			size_t ns = lrs_.left().permutationVector().size();
			//size_t nx = ns/A.data.rank();
			if (src.size()!=lrs_.super().permutationVector().size())
				throw std::runtime_error("applyLocalOpSystem SE\n");

			PackIndicesType pack(ns);

			for (size_t i=offset;i<final;i++) {
				size_t x=0,y=0;
				pack.unpack(x,y,lrs_.super().permutation(i));
				//if (y>=basisE_.permutationVector().size()) throw std::runtime_error("applyLocalOpSystem E\n");
				if (x>=lrs_.left().permutationVector().size())
					throw std::runtime_error("applyLocalOpSystem S\n");
				RealType sign = lrs_.left().fermionicSign(x,A.fermionSign);
				for (int k=A.data.getRowPtr(y);k<A.data.getRowPtr(y+1);k++) {
					size_t yprime = A.data.getCol(k);
					size_t j = lrs_.super().permutationInverse(x+yprime*ns);
					dest2[j] += src[i]*A.data.getValue(k)*sign;
				}
			}
		}

		// entry point for corner cases. These are all when expanding ths system
		void applyLocalOpCorner(
				VectorWithOffsetType& dest,
				const VectorWithOffsetType& src,
				const OperatorType& A,
				const FermionSign& fermionSign) const
		{
			if (lrs_.right().size() == A.data.rank()) { // right corner
				applyLocalOpSystem(dest,src,A,fermionSign,RIGHT_CORNER);
				return;
			}
			applyLocalOpEnviron(dest,src,A,LEFT_CORNER);
		}

		const LeftRightSuperType& lrs_;
	}; // class ApplyOperatorLocal
} // namespace Dmrg 

/*@}*/
#endif //APPLY_OPERATOR_LOCAL_H
