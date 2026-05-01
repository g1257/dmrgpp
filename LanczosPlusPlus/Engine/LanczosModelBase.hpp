
/*
// BEGIN LICENSE BLOCK
Copyright (c) 2014, UT-Battelle, LLC
All rights reserved

[Lanczos++, Version 1.0.0]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

*/

#ifndef LANCZOS_MODEL_BASE_H
#define LANCZOS_MODEL_BASE_H
#include "BasisBase.h"
#include "CrsMatrix.h"
#include "RahulOperator.h"
#include "Vector.h"

namespace LanczosPlusPlus {

template <typename ComplexOrRealType_, typename GeometryType_, typename InputType_>

class LanczosModelBase {

public:

	typedef typename PsimagLite::Real<ComplexOrRealType_>::Type  RealType;
	typedef GeometryType_                                        GeometryType;
	typedef ComplexOrRealType_                                   ComplexOrRealType;
	typedef InputType_                                           InputType;
	typedef BasisBase<GeometryType>                              BasisBaseType;
	typedef PsimagLite::CrsMatrix<ComplexOrRealType>             SparseMatrixType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef typename PsimagLite::Vector<SizeType>::Type          VectorSizeType;
	typedef RahulOperator<ComplexOrRealType>                     RahulOperatorType;
	typedef typename PsimagLite::Vector<RahulOperatorType>::Type VectorRahulOperatorType;
	typedef typename BasisBaseType::WordType                     WordType;

	virtual ~LanczosModelBase() { }

	virtual SizeType size() const = 0;

	virtual SizeType orbitals(SizeType site) const = 0;

	virtual void setupHamiltonian(SparseMatrixType& matrix) const = 0;

	virtual bool hasNewParts(std::pair<SizeType, SizeType>&       newParts,
	                         const std::pair<SizeType, SizeType>& oldParts,
	                         const LabeledOperator&               what,
	                         SizeType                             spin,
	                         SizeType                             orb) const
	    = 0;

	virtual const GeometryType& geometry() const = 0;

	virtual void setupHamiltonian(SparseMatrixType& matrix, const BasisBaseType& basis) const
	    = 0;

	virtual void matrixVectorProduct(VectorType&, const VectorType&) const
	{
		throw PsimagLite::RuntimeError(
		    "ModelBase::matrixVectorProduct(2) not impl. for this model\n");
	}

	virtual void matrixVectorProduct(VectorType&, const VectorType&, const BasisBaseType&) const
	{
		throw PsimagLite::RuntimeError(
		    "ModelBase::matrixVectorProduct(3) not impl. for this model\n");
	}

	virtual const BasisBaseType& basis() const = 0;

	virtual PsimagLite::String name() const = 0;

	virtual BasisBaseType* createBasis(SizeType nup, SizeType ndown) const = 0;

	virtual void print(std::ostream& os) const = 0;

	void rahulMethod(VectorType&                    psiNew,
	                 const VectorRahulOperatorType& vops,
	                 const VectorSizeType&          vsites,
	                 const VectorType&              psi,
	                 const BasisBaseType&           basis) const
	{
		static const int FERMION_SIGN = LanczosGlobals::FERMION_SIGN;

		std::fill(psiNew.begin(), psiNew.end(), 0.0);
		const SizeType hilbert = basis.size();
		const SizeType nops    = vops.size();
		for (SizeType ispace = 0; ispace < hilbert; ++ispace) {
			WordType ket1 = basis(ispace, LanczosGlobals::SPIN_UP);
			WordType ket2 = basis(ispace, LanczosGlobals::SPIN_DOWN);
			assert(ispace < psi.size());
			ComplexOrRealType value   = psi[ispace];
			WordType          ketp1   = ket1;
			WordType          ketp2   = ket2;
			bool              nonZero = false;
			for (SizeType jj = 0; jj < nops; ++jj) {
				const SizeType           j  = nops - jj - 1; // start from the end
				const RahulOperatorType& op = vops[j];
				assert(j < vsites.size());
				const SizeType site = vsites[j];

				ComplexOrRealType result = 0;
				nonZero                  = (op.dof() == LanczosGlobals::SPIN_UP)
				                     ? applyOperator(ketp1, result, op, site)
				                     : applyOperator(ketp2, result, op, site);
				if (!nonZero)
					break;

				if (op.isFermionic()) {
					const SizeType overUp
					    = (op.dof() == LanczosGlobals::SPIN_UP)
					    ? 0
					    : PsimagLite::BitManip::count(ketp1);

					if (overUp & 1)
						result *= FERMION_SIGN;

					const WordType ket1or2
					    = (op.dof() == LanczosGlobals::SPIN_UP) ? ketp1 : ketp2;
					const ComplexOrRealType fsign
					    = LanczosGlobals::doSign(ket1or2, site);
					result *= fsign;
				}

				value *= result;
			}

			if (!nonZero)
				continue;
			SizeType newI = basis.perfectIndex(ketp1, ketp2);
			assert(newI < psiNew.size());
			psiNew[newI] += value;
		}
	}

	virtual void printOperators(std::ostream&) const
	{
		PsimagLite::String str(__FILE__);
		str += " " + ttos(__LINE__) + "\n";
		str += PsimagLite::String("Function printOperators unimplemented\n");
		std::cerr << str.c_str();
	}

protected:

	template <typename SomeVectorType>
	static typename PsimagLite::EnableIf<PsimagLite::IsVectorLike<SomeVectorType>::True,
	                                     void>::Type
	deleteGarbage(SomeVectorType& garbage)
	{
		for (SizeType i = 0; i < garbage.size(); ++i) {
			delete garbage[i];
			garbage[i] = 0;
		}
	}

private:

	static bool applyOperator(WordType&                ketp,
	                          ComplexOrRealType&       result,
	                          const RahulOperatorType& op,
	                          SizeType                 site)
	{
		const WordType ket       = ketp;
		const WordType mask      = LanczosGlobals::bitmask(site);
		WordType       s         = (ket & mask);
		bool           sbit      = (s > 0);
		bool           sbitSaved = sbit;
		bool           nonZero   = op.actOn(sbit, result);
		if (!nonZero)
			return false;
		if (sbitSaved != sbit)
			ketp ^= mask;
		return true;
	}

}; // class ModelBase

template <typename RealType, typename GeometryType, typename InputType>
std::ostream& operator<<(std::ostream&                                              os,
                         const LanczosModelBase<RealType, GeometryType, InputType>& model)
{
	model.print(os);
	return os;
}

} // namespace LanczosPlusPlus

#endif
