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

/** \ingroup DMRG */
/*@{*/

/*! \file ApplyOperatorLocal.h
 *
 *  dest2 = transpose(A) * src
 */
#ifndef APPLY_OPERATOR_LOCAL_H
#define APPLY_OPERATOR_LOCAL_H

#include "PackIndices.h" // in PsimagLite
#include "FermionSign.h"
#include "ProgramGlobals.h"
#include "OperatorStorage.h"

namespace Dmrg {

template<typename LeftRightSuperType_, typename VectorWithOffsetType_>
class ApplyOperatorLocal {

	typedef typename VectorWithOffsetType_::VectorType TargetVectorType;
	typedef typename LeftRightSuperType_::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisWithOperatorsType::ComplexOrRealType ComplexOrRealType;
	typedef PsimagLite::PackIndices PackIndicesType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType_;

	class LegacyBug {

	public:

		LegacyBug(bool withLegacyBug, const OperatorType_& A)
		    : withLegacyBug_(withLegacyBug),
		      Aptr_((withLegacyBug) ? const_cast<OperatorType_*>(&A ):
		                              new OperatorType_)
		{
			if (withLegacyBug_) return;

			*Aptr_ = A;

			Aptr_->transpose();
		}

		~LegacyBug()
		{
			if (withLegacyBug_) return;
			delete Aptr_;
			Aptr_ = 0;
		}

		const OperatorType_& operator()() const
		{
			assert(Aptr_);
			return *Aptr_;
		}

	private:

		bool withLegacyBug_;
		OperatorType_* Aptr_;
	}; // class LegacyBug

public:

	enum BorderEnum {BORDER_NO = false, BORDER_YES = true};

	enum {MIDDLE, LEFT_CORNER, RIGHT_CORNER};

	typedef LeftRightSuperType_ LeftRightSuperType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef VectorWithOffsetType_ VectorWithOffsetType;
	typedef OperatorType_ OperatorType;
	typedef FermionSign FermionSignType;

	ApplyOperatorLocal(const LeftRightSuperType& lrs, bool withLegacyBug)
	    : lrs_(lrs), withLegacyBug_(withLegacyBug)
	{}

	//! FIXME: we need to make a fast version for when we're just
	//! figuring out where the (non-zero) partition is
	void operator()(VectorWithOffsetType& dest,
	                const VectorWithOffsetType& src,
	                const OperatorType& AA,
	                const FermionSign& fermionSign,
	                ProgramGlobals::DirectionEnum systemOrEnviron,
	                BorderEnum corner) const
	{
		LegacyBug legacyBug(withLegacyBug_, AA);
		const OperatorType& A = legacyBug();

		if (corner == BORDER_NO) {
			if (systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
				applyLocalOpSystem(dest,src,A,fermionSign);
			else
				applyLocalOpEnviron(dest,src,A);
			return;
		}

		applyLocalOpCorner(dest,src,A,fermionSign);
	}

	//! FIXME: we need to make a fast version for when we're just
	//! figuring out where the (non-zero) partition is
	void hookForZero(VectorWithOffsetType& dest,
	                 const VectorWithOffsetType& src,
	                 const OperatorType& A,
	                 const FermionSign& fermionSign,
	                 ProgramGlobals::DirectionEnum systemOrEnviron) const
	{
		assert(systemOrEnviron == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM);

		TargetVectorType dest2(lrs_.super().size(),0.0);

		for (SizeType ii = 0; ii < src.sectors(); ++ii) {
			SizeType i = src.sector(ii);
			hookForZeroSystem(dest2,src,A,fermionSign,i);
		}

		dest.fromFull(dest2,lrs_.super());
	}

	// dest2 = transpose(A) * src; corrected if !withLegacyBug
	void hookForZeroSystem(TargetVectorType& dest2,
	                       const VectorWithOffsetType& src,
	                       const OperatorType& AA,
	                       const FermionSign&,
	                       SizeType i0) const
	{
		LegacyBug legacyBug(withLegacyBug_, AA);
		const OperatorType& A = legacyBug();

		SizeType offset = src.offset(i0);
		SizeType final = offset + src.effectiveSize(i0);
		SizeType ns = lrs_.left().permutationVector().size();
		SizeType nx = ns/A.getCRS().rows();
		if (src.size() != lrs_.super().permutationVector().size())
			err("applyLocalOpSystem SE\n");

		PackIndicesType pack1(ns);
		PackIndicesType pack2(nx);
		for (SizeType i = offset; i < final; ++i) {
			SizeType x = 0;
			SizeType y  = 0;
			pack1.unpack(x,y,lrs_.super().permutation(i));

			SizeType x0 = 0;
			SizeType x1 = 0;
			assert(x<lrs_.left().permutationVector().size());
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			SizeType start = A.getCRS().getRowPtr(x0);
			SizeType end = A.getCRS().getRowPtr(x0 + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType x0prime = A.getCRS().getCol(k);
				SizeType xprime = lrs_.left().permutationInverse(x0prime+x1*nx);
				SizeType j = lrs_.super().permutationInverse(xprime+y*ns);
				dest2[j] += src.slowAccess(i)*A.getCRS().getValue(k);
			}
		}
	}

	const LeftRightSuperType& lrs() const { return lrs_; }

private:

	ApplyOperatorLocal(const ApplyOperatorLocal&);

	ApplyOperatorLocal& operator=(const ApplyOperatorLocal&);

	void applyLocalOpSystem(VectorWithOffsetType& dest,
	                        const VectorWithOffsetType& src,
	                        const OperatorType& A,
	                        const FermionSign& fermionSign,
	                        SizeType whichPartOfTheLattice = MIDDLE) const
	{
		TargetVectorType dest2(lrs_.super().size(),0.0);

		for (SizeType ii=0;ii<src.sectors();ii++) {
			SizeType i = src.sector(ii);
			switch (whichPartOfTheLattice) {
			case MIDDLE:
				applyLocalOpSystem(dest2,src,A,fermionSign,i);
				break;
			case LEFT_CORNER:
				throw PsimagLite::RuntimeError("applyLocalOpSystem: internal error\n");
			case RIGHT_CORNER:
				applyLocalOpRightCorner(dest2,src,A,i);
				break;
			}
		}

		dest.fromFull(dest2,lrs_.super());
	}

	// dest2 = transpose(A) * src; corrected if !withLegacyBug
	void applyLocalOpSystem(TargetVectorType& dest2,
	                        const VectorWithOffsetType& src,
	                        const OperatorType& A,
	                        const FermionSign& fermionSign,
	                        SizeType i0) const
	{
		SizeType offset = src.offset(i0);
		SizeType final = offset + src.effectiveSize(i0);
		SizeType ns = lrs_.left().permutationVector().size();
		SizeType nx = ns/A.getCRS().rows();
		if (src.size() != lrs_.super().permutationVector().size())
			err("applyLocalOpSystem SE\n");

		PackIndicesType pack1(ns);
		PackIndicesType pack2(nx);
		for (SizeType i = offset; i < final; ++i) {
			SizeType x = 0;
			SizeType y = 0;
			pack1.unpack(x,y,lrs_.super().permutation(i));

			SizeType x0 = 0;
			SizeType x1 = 0;
			assert(x<lrs_.left().permutationVector().size());
			pack2.unpack(x0,x1,lrs_.left().permutation(x));

			const bool isFermion = (A.fermionOrBoson() ==
			                        ProgramGlobals::FermionOrBosonEnum::FERMION);
			const RealType sign = fermionSign(x0, (isFermion) ? -1 : 1);
			const SizeType start = A.getCRS().getRowPtr(x1);
			const SizeType end = A.getCRS().getRowPtr(x1 + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType x1prime = A.getCRS().getCol(k);
				SizeType xprime = lrs_.left().permutationInverse(x0+x1prime*nx);
				SizeType j = lrs_.super().permutationInverse(xprime+y*ns);
				dest2[j] += src.slowAccess(i)*A.getCRS().getValue(k)*sign;
			}
		}
	}

	void applyLocalOpEnviron(VectorWithOffsetType& dest,
	                         const VectorWithOffsetType& src,
	                         const OperatorType& A,
	                         SizeType whichPartOfTheLattice = MIDDLE) const
	{
		TargetVectorType dest2(lrs_.super().size(),0.0);

		for (SizeType ii=0;ii<src.sectors();ii++) {
			SizeType i = src.sector(ii);
			switch (whichPartOfTheLattice) {
			case MIDDLE:
				applyLocalOpEnviron(dest2,src,A,i);
				break;
			case LEFT_CORNER:
				applyLocalOpLeftCorner(dest2,src,A,i);
				break;
			case RIGHT_CORNER:
				err("applyLocalOpEnviron: internal error\n");
			}
		}

		dest.fromFull(dest2,lrs_.super());
	}

	// dest2 = transpose(A) * src; corrected if !withLegacyBug
	void applyLocalOpEnviron(TargetVectorType& dest2,
	                         const VectorWithOffsetType& src,
	                         const OperatorType& A,
	                         SizeType i0) const
	{
		SizeType offset = src.offset(i0);
		SizeType final = offset + src.effectiveSize(i0);

		SizeType ns = lrs_.left().size();
		SizeType nx = A.getCRS().rows();
		PackIndicesType pack1(ns);
		PackIndicesType pack2(nx);

		for (SizeType i = offset; i < final; ++i) {
			SizeType x = 0;
			SizeType y = 0;
			pack1.unpack(x,y,lrs_.super().permutation(i));
			SizeType y0 = 0;
			SizeType y1 = 0;
			pack2.unpack(y0,y1,lrs_.right().permutation(y));
			const bool isFermion = (A.fermionOrBoson() ==
			                        ProgramGlobals::FermionOrBosonEnum::FERMION);
			const RealType sign = lrs_.left().fermionicSign(x, (isFermion) ? -1 : 1);
			const SizeType start = A.getCRS().getRowPtr(y0);
			const SizeType end = A.getCRS().getRowPtr(y0 + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType y0prime = A.getCRS().getCol(k);
				SizeType yprime = lrs_.right().permutationInverse(y0prime+y1*nx);
				SizeType j = lrs_.super().permutationInverse(x+yprime*ns);
				dest2[j] += src.slowAccess(i)*A.getCRS().getValue(k)*sign;
			}
		}
	}

	// dest2 = transpose(A) * src; corrected if !withLegacyBug
	void applyLocalOpLeftCorner(TargetVectorType& dest2,
	                            const VectorWithOffsetType& src,
	                            const OperatorType& A,
	                            SizeType i0) const
	{
		SizeType offset = src.offset(i0);
		SizeType final = offset + src.effectiveSize(i0);

		SizeType ns = lrs_.left().size();
		PackIndicesType pack(ns);

		for (SizeType i = offset; i < final; ++i) {
			SizeType x = 0;
			SizeType y = 0;
			pack.unpack(x,y,lrs_.super().permutation(i));

			const SizeType start = A.getCRS().getRowPtr(x);
			const SizeType end = A.getCRS().getRowPtr(x + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType xprime = A.getCRS().getCol(k);
				SizeType j = lrs_.super().permutationInverse(xprime+y*ns);
				dest2[j] += src.slowAccess(i)*A.getCRS().getValue(k);
			}
		}
	}

	// dest2 = transpose(A) * src; corrected if !withLegacyBug
	void applyLocalOpRightCorner(TargetVectorType& dest2,
	                             const VectorWithOffsetType& src,
	                             const OperatorType& A,
	                             SizeType i0) const
	{
		SizeType offset = src.offset(i0);
		SizeType final = offset + src.effectiveSize(i0);
		SizeType ns = lrs_.left().permutationVector().size();
		if (src.size() != lrs_.super().permutationVector().size())
			err("applyLocalOpSystem SE\n");

		PackIndicesType pack(ns);

		for (SizeType i = offset; i < final; ++i) {
			SizeType x = 0;
			SizeType y = 0;
			pack.unpack(x,y,lrs_.super().permutation(i));

			if (x >= lrs_.left().permutationVector().size())
				err("applyLocalOpSystem S\n");

			const bool isFermion = (A.fermionOrBoson() ==
			                        ProgramGlobals::FermionOrBosonEnum::FERMION);
			const RealType sign = lrs_.left().fermionicSign(x, (isFermion) ? -1 : 1);
			const SizeType start = A.getCRS().getRowPtr(y);
			const SizeType end = A.getCRS().getRowPtr(y + 1);
			for (SizeType k = start; k < end; ++k) {
				SizeType yprime = A.getCRS().getCol(k);
				SizeType j = lrs_.super().permutationInverse(x+yprime*ns);
				dest2[j] += src.slowAccess(i)*A.getCRS().getValue(k)*sign;
			}
		}
	}

	// entry point for corner cases. These are all when expanding ths system
	void applyLocalOpCorner(VectorWithOffsetType& dest,
	                        const VectorWithOffsetType& src,
	                        const OperatorType& A,
	                        const FermionSign& fermionSign) const
	{
		if (lrs_.right().size() == A.getCRS().rows()) { // right corner
			applyLocalOpSystem(dest,src,A,fermionSign,RIGHT_CORNER);
			return;
		}

		applyLocalOpEnviron(dest,src,A,LEFT_CORNER);
	}

	const LeftRightSuperType& lrs_;
	bool withLegacyBug_;
}; // class ApplyOperatorLocal
} // namespace Dmrg

/*@}*/
#endif //APPLY_OPERATOR_LOCAL_H

