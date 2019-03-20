/*
Copyright (c)  2008-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 5.]
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

/*! \file CorrelationsSkeleton.h
 *
 *  Helper class for Observables
 *
 */
#ifndef CORRELATIONS_SK_H
#define CORRELATIONS_SK_H
#include "ObserverHelper.h"
#include "Matrix.h"
#include "PackIndices.h"
#include "CrsMatrix.h"
#include "ApplyOperatorLocal.h"
#include "Braket.h"
#include <numeric>

namespace Dmrg {

template<typename ObserverHelperType_,typename ModelType>
class CorrelationsSkeleton {

public:

	typedef PsimagLite::PackIndices PackIndicesType;
	typedef Braket<ModelType> BraketType;
	typedef ObserverHelperType_ ObserverHelperType;
	typedef typename ObserverHelperType::MatrixType MatrixType;
	typedef typename ObserverHelperType::VectorType VectorType ;
	typedef typename ObserverHelperType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType BasisWithOperatorsType ;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename ObserverHelperType::FermionSignType FermionSignType;
	typedef typename BasisType::VectorSizeType VectorSizeType;
	typedef typename VectorType::value_type FieldType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef PsimagLite::CrsMatrix<FieldType> SparseMatrixType;

	enum class GrowDirection {RIGHT, LEFT};

	CorrelationsSkeleton(const ObserverHelperType& helper) : helper_(helper)
	{}

	SizeType numberOfSites() const
	{
		return helper_.numberOfSites();
	}

	//! i can be zero here!!
	void growDirectly(SparseMatrixType& Odest,
	                  const SparseMatrixType& Osrc,
	                  SizeType i,
	                  ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                  SizeType ns,
	                  bool transform) const
	{
		Odest =Osrc;
		// from 0 --> i
		int nt=i-1;
		if (nt<0) nt=0;

		for (SizeType s = nt; s < ns; ++s) {
			const GrowDirection growOption = growthDirection(s, nt, i, s);
			SparseMatrixType Onew(helper_.cols(s),helper_.cols(s));

			fluffUp(Onew, Odest, fermionicSign, growOption, false, s);
			if (!transform && s == ns - 1) {
				Odest = Onew;
				continue;
			}

			helper_.transform(Odest, Onew, s);
		}
	}

	GrowDirection growthDirection(SizeType s,
	                              int nt,
	                              SizeType i,
	                              SizeType ptr) const
	{
		const ProgramGlobals::DirectionEnum dir = helper_.direction(ptr);
		GrowDirection growOption = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? GrowDirection::RIGHT
		        : GrowDirection::LEFT;

		if (s==SizeType(nt)) {
			growOption = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			        ? GrowDirection::LEFT
			        : GrowDirection::RIGHT;
			if (i==0) growOption = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			        ? GrowDirection::RIGHT
			        : GrowDirection::LEFT;
		}

		return growOption;
	}

	// Perfomance critical:
	void fluffUp(SparseMatrixType& ret2,
	             const SparseMatrixType& O,
	             ProgramGlobals::FermionOrBosonEnum fOrB,
	             const GrowDirection growOption,
	             bool transform,
	             SizeType ptr) const
	{
		const int fermionicSign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;
		const ProgramGlobals::DirectionEnum dir = helper_.direction(ptr);

		const BasisType& basis = (dir == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? helper_.leftRightSuper(ptr).left()
		        : helper_.leftRightSuper(ptr).right();

		SizeType n = basis.size();
		SizeType orows = O.rows();
		SizeType ktotal = n/orows;
		SparseMatrixType ret(n, n, ktotal*O.nonZeros());

		if (growOption == GrowDirection::RIGHT) {
			RealType sign = (dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON)
			        ? fermionSignBasis(fermionicSign, helper_.leftRightSuper(ptr).left()) :
			          1;

			SizeType counter = 0;
			for (SizeType k = 0; k < ktotal; ++k) {
				for (SizeType i = 0; i < orows; ++i) {
					ret.setRow(i + k*orows, counter);
					for (int kj = O.getRowPtr(i); kj < O.getRowPtr(i + 1); ++kj) {
						// Sperm[e0] = i + k*n
						// Sperm[e1] = j + k*n

						SizeType col = basis.
						        permutationInverse(O.getCol(kj) + k*orows);

						ret.setCol(counter, col);
						ret.setValues(counter++, O.getValue(kj)*sign);
					}
				}
			}

			ret.setRow(n, counter);
			ret.checkValidity();
		} else {

			SizeType counter = 0;
			for (SizeType i = 0; i < orows; ++i) {
				for (SizeType k = 0; k < ktotal; ++k) {
					ret.setRow(k + i*ktotal, counter);
					RealType sign = 1;
					if (dir == ProgramGlobals::DirectionEnum::EXPAND_ENVIRON) {
						sign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON)
						        ? 1 : fermionSignBasis(fermionicSign,
						                               helper_.leftRightSuper(ptr).left())*
						          helper_.signsOneSite(k);
					} else {
						sign = helper_.fermionicSignLeft(ptr)(k, fermionicSign);
					}

					for (int kj = O.getRowPtr(i); kj < O.getRowPtr(i + 1); ++kj) {
						// Sperm[e0] = k + i*m
						// Sperm[e1] = k + j*m

						SizeType col = basis.
						        permutationInverse(k + O.getCol(kj)*ktotal);

						ret.setCol(counter, col);
						ret.setValues(counter++, O.getValue(kj)*sign);
					}
				}
			}

			ret.setRow(n, counter);
			ret.checkValidity();
		}

		reorderRowsCrs(ret, basis);
		if (transform) {
			helper_.transform(ret2, ret, ptr);
			return;
		}

		ret2 = ret;
	}

	SizeType dmrgMultiply(SparseMatrixType& result,
	                      const SparseMatrixType& O1,
	                      const SparseMatrixType& O2,
	                      ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                      SizeType ns) const
	{
		SizeType ptr = (ns == 0) ? ns : ns - 1;
		if (helper_.direction(ptr) == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			dmrgMultiplySystem(result, ptr, O1, O2, fermionicSign, ns);
		else
			dmrgMultiplyEnviron(result, ptr, O1, O2, fermionicSign, ns);

		return ptr;

	}

	static void createWithModification(SparseMatrixType& Om,const SparseMatrixType& O,char mod)
	{
		if (mod == 'n' || mod == 'N') {
			Om = O;
			return;
		}

		transposeConjugate(Om,O);
	}

	FieldType bracket(const SparseMatrixType& A,
	                  ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                  SizeType ptr,
	                  PsimagLite::String bra,
	                  PsimagLite::String ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);

			return bracket_(A,src1,src2,fermionicSign,ptr);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: CorrelationsSkeleton::bracket(...):";
			std::cerr<<" No data seen yet\n";
			return 0;
		}
	}

	FieldType bracketRightCorner(const SparseMatrixType& A,
	                             const SparseMatrixType& B,
	                             ProgramGlobals::FermionOrBosonEnum fermionSign,
	                             SizeType ptr,
	                             PsimagLite::String bra,
	                             PsimagLite::String ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);
			return bracketRightCorner_(A,B,fermionSign,src1,src2,ptr);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: CorrelationsSkeleton::bracketRightCorner(...):";
			std::cerr<<" No data seen yet\n";
			return 0;
		}
	}

	FieldType bracketRightCorner(const SparseMatrixType& A,
	                             const SparseMatrixType& B,
	                             const SparseMatrixType& C,
	                             ProgramGlobals::FermionOrBosonEnum fermionSign,
	                             SizeType ptr,
	                             PsimagLite::String bra,
	                             PsimagLite::String ket) const
	{
		try {
			const VectorWithOffsetType& src1 = helper_.getVectorFromBracketId(bra, ptr);
			const VectorWithOffsetType& src2 = helper_.getVectorFromBracketId(ket, ptr);
			return bracketRightCorner_(A,B,C,fermionSign,src1,src2,ptr);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: CorrelationsSkeleton::bracketRightCornerABC(...):";
			std::cerr<<" No data seen yet\n";
			return 0;
		}
	}

	const ObserverHelperType& helper() const { return helper_; }

private:

	int fermionSignBasis(int fermionicSign, const BasisType& basis) const
	{
		const typename BasisWithOperatorsType::VectorBoolType& v = basis.signs();
		SizeType n = v.size();
		if (n == 0) return 1;
		bool nx0 = v[0];
		for (SizeType i = 1; i < n; ++i)
			nx0 ^= v[i];
		return (nx0) ? fermionicSign : 1;
	}

	void dmrgMultiplySystem(SparseMatrixType& result,
	                        SizeType& ptr,
	                        const SparseMatrixType& O1,
	                        const SparseMatrixType& O2,
	                        ProgramGlobals::FermionOrBosonEnum fOrB, // for O2
	                        SizeType ns) const
	{
		const int fermionicSign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;
		SizeType ni=O1.rows();

		ptr = ns;
		SizeType sprime = helper_.leftRightSuper(ptr).left().size(); //ni*nj;
		result.resize(sprime,sprime);

		if (helper_.leftRightSuper(ptr).left().size()!=sprime) {
			std::cerr<<"WARNING: "<<helper_.leftRightSuper(ptr).left().size();
			std::cerr<<"!="<<sprime<<"\n";
			err("problem in dmrgMultiply\n");
		}

		PsimagLite::Vector<SizeType>::Type col(sprime,0);
		typename PsimagLite::Vector<FieldType>::Type value(sprime, 0);

		PackIndicesType pack(ni);

		SizeType counter = 0;
		for (SizeType r=0;r<sprime;r++) {
			SizeType e,u;
			pack.unpack(e,u,helper_.leftRightSuper(ptr).left().permutation(r));
			RealType f = helper_.fermionicSignLeft(ptr)(e,fermionicSign);
			result.setRow(r,counter);
			for (int k=O1.getRowPtr(e);k<O1.getRowPtr(e+1);k++) {
				SizeType e2 = O1.getCol(k);
				for (int k2=O2.getRowPtr(u);k2<O2.getRowPtr(u+1);k2++) {
					SizeType u2 = O2.getCol(k2);
					SizeType r2 = helper_.leftRightSuper(ptr).left().
					        permutationInverse(e2 + u2*ni);
					value[r2] += O1.getValue(k)*O2.getValue(k2)*f;
					col[r2] = 1;
				}
			}

			for (SizeType i=0;i<col.size();i++) {
				if (col[i]==0) continue;
				result.pushCol(i);
				result.pushValue(value[i]);
				counter++;
				value[i] = 0.0;
				col[i] = 0;
			}
		}

		result.setRow(result.rows(),counter);
		result.checkValidity();
	}

	void dmrgMultiplyEnviron(SparseMatrixType& result,
	                         SizeType& ptr,
	                         const SparseMatrixType& O1,
	                         const SparseMatrixType& O2,
	                         ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                         SizeType ns) const
	{
		int fs = (fermionicSign == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;
		RealType f = fermionSignBasis(fs,
		                              helper_.leftRightSuper(ptr).left());
		SizeType nj=O2.rows();

		ptr = ns;
		SizeType eprime = helper_.leftRightSuper(ptr).right().size(); //ni*nj;
		result.resize(eprime,eprime);

		if (helper_.leftRightSuper(ptr).right().size()!=eprime) {
			std::cerr<<"WARNING: "<<helper_.leftRightSuper(ptr).right().size();
			std::cerr<<"!="<<eprime<<"\n";
			throw PsimagLite::RuntimeError("problem in dmrgMultiply\n");
		}

		PsimagLite::Vector<SizeType>::Type col(eprime,0);
		typename PsimagLite::Vector<FieldType>::Type value(eprime,0);

		PackIndicesType pack(nj);

		SizeType counter = 0;
		for (SizeType r=0;r<eprime;r++) {
			result.setRow(r,counter);
			SizeType e = 0;
			SizeType u = 0;
			pack.unpack(e,u,helper_.leftRightSuper(ptr).right().permutation(r));
			const RealType sign = (fermionicSign == ProgramGlobals::FermionOrBosonEnum::BOSON)
			        ? 1 : f*helper_.signsOneSite(e);

			for (int k=O2.getRowPtr(e);k<O2.getRowPtr(e+1);k++) {
				SizeType e2 = O2.getCol(k);

				for (int k2=O1.getRowPtr(u);k2<O1.getRowPtr(u+1);k2++) {
					SizeType u2 = O1.getCol(k2);
					SizeType r2 = helper_.leftRightSuper(ptr).right().
					        permutationInverse(e2 + u2*nj);
					assert(r2<eprime);
					col[r2] = 1;
					value[r2] += O2.getValue(k)*O1.getValue(k2)*sign;
				}
			}

			for (SizeType i=0;i<col.size();i++) {
				if (col[i]==0) continue;
				result.pushCol(i);
				result.pushValue(value[i]);
				counter++;
				value[i] = 0.0;
				col[i] = 0;
			}
		}

		result.setRow(result.rows(),counter);
		result.checkValidity();
	}

	static void reorderRowsCrs(SparseMatrixType& matrix, const BasisType& basis)
	{
		SizeType nrows = matrix.rows();
		SizeType ncols = matrix.cols();
		SparseMatrixType matrix2(nrows, ncols, matrix.nonZeros());

		SizeType counter = 0;
		for (SizeType row = 0; row < nrows; ++row) {
			matrix2.setRow(row, counter);
			SizeType r = basis.permutation(row);
			SizeType start = matrix.getRowPtr(r);
			SizeType end = matrix.getRowPtr(r + 1);
			for (SizeType k = start; k < end; ++k) {
				matrix2.setCol(counter, matrix.getCol(k));
				matrix2.setValues(counter++, matrix.getValue(k));
			}
		}

		matrix2.setRow(nrows, counter);
		matrix2.checkValidity();
		matrix = matrix2;
	}

	FieldType bracket_(const SparseMatrixType& A,
	                   const VectorWithOffsetType& vec1,
	                   const VectorWithOffsetType& vec2,
	                   ProgramGlobals::FermionOrBosonEnum fermionicSign,
	                   SizeType ptr) const
	{
		if (vec1.size()!=helper_.leftRightSuper(ptr).super().size() ||
		        vec1.size()!=vec2.size())
			err("CorrelationsSkeleton::bracket_(...): Error\n");

		return (helper_.direction(ptr) == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? bracketSystem_(A,vec1,vec2,ptr)
		        : bracketEnviron_(A,vec1,vec2,fermionicSign,ptr);
	}

	FieldType bracketSystem_(const SparseMatrixType& A,
	                         const VectorWithOffsetType& vec1,
	                         const VectorWithOffsetType& vec2,
	                         SizeType ptr) const
	{
		FieldType sum=0;
		PackIndicesType pack(helper_.leftRightSuper(ptr).left().size());
		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack.unpack(r,eta,helper_.leftRightSuper(ptr).super().
				            permutation(t));
				for (int k=A.getRowPtr(r);k<A.getRowPtr(r+1);k++) {
					SizeType r2 = A.getCol(k);
					SizeType t2 = helper_.leftRightSuper(ptr).super().
					        permutationInverse(r2+eta*A.cols());
					if (t2<offset || t2>=total) continue;
					sum += A.getValue(k)*PsimagLite::conj(vec1.slowAccess(t))*
					        vec2.slowAccess(t2);
				}
			}
		}

		return resultDivided(sum,vec1);
	}

	FieldType bracketEnviron_(const SparseMatrixType& A,
	                          const VectorWithOffsetType& vec1,
	                          const VectorWithOffsetType& vec2,
	                          ProgramGlobals::FermionOrBosonEnum fOrB,
	                          SizeType ptr) const
	{
		const int fermionicSign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;

		RealType sign = fermionSignBasis(fermionicSign,
		                                 helper_.leftRightSuper(ptr).left());

		FieldType sum=0;
		PackIndicesType pack(helper_.leftRightSuper(ptr).left().size());
		SizeType leftSize = helper_.leftRightSuper(ptr).left().size();

		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack.unpack(r,eta,helper_.leftRightSuper(ptr).super().
				            permutation(t));
				if (eta>=A.rows()) throw PsimagLite::RuntimeError("Error\n");

				for (int k=A.getRowPtr(eta);k<A.getRowPtr(eta+1);k++) {
					SizeType eta2 = A.getCol(k);
					SizeType t2 = helper_.leftRightSuper(ptr).super().
					        permutationInverse(r+eta2*leftSize);
					if (t2<offset || t2>=total) continue;
					sum += A.getValue(k)*PsimagLite::conj(vec1.slowAccess(t))*
					        vec2.slowAccess(t2)*sign;
				}
			}
		}

		return resultDivided(sum,vec1);
	}

	FieldType bracketRightCorner_(const SparseMatrixType& A,
	                              const SparseMatrixType& B,
	                              ProgramGlobals::FermionOrBosonEnum fermionSign,
	                              const VectorWithOffsetType& vec1,
	                              const VectorWithOffsetType& vec2,
	                              SizeType ptr) const
	{
		return (helper_.direction(ptr) == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
		        ? brRghtCrnrSystem_(A,B,fermionSign,vec1,vec2,ptr)
		        : brLftCrnrEnviron_(A,B,fermionSign,vec1,vec2,ptr);
	}

	bool superOddElectrons(SizeType t, SizeType ptr) const
	{
		// return helper_.leftRightSuper(ptr).super().electrons(t);
		SizeType tmp = helper_.leftRightSuper(ptr).super().permutation(t);
		div_t mydiv = PsimagLite::div(tmp,helper_.leftRightSuper(ptr).left().size());
		return helper_.leftRightSuper(ptr).right().signs()[mydiv.quot] ^
		        helper_.leftRightSuper(ptr).left().signs()[mydiv.rem];
	}

	FieldType brRghtCrnrSystem_(const SparseMatrixType& Acrs,
	                            const SparseMatrixType& Bcrs,
	                            ProgramGlobals::FermionOrBosonEnum fermionSign,
	                            const VectorWithOffsetType& vec1,
	                            const VectorWithOffsetType& vec2,
	                            SizeType ptr) const
	{
		FieldType sum = 0;
		SizeType leftSize = helper_.leftRightSuper(ptr).left().size();
		SizeType ni = helper_.leftRightSuper(ptr).left().size()/Bcrs.rows();

		// some sanity checks:
		if (vec1.size() != vec2.size() ||
		        vec1.size() != helper_.leftRightSuper(ptr).super().size())
			err("Observe::brRghtCrnrSystem_(...)\n");

		if (ni != Acrs.rows())
			err("Observe::brRghtCrnrSystem_(...)\n");

		// ok, we're ready for the main course:
		PackIndicesType pack1(helper_.leftRightSuper(ptr).left().size());
		PackIndicesType pack2(ni);
		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack1.unpack(r,eta,helper_.leftRightSuper(ptr).super().
				             permutation(t));
				SizeType r0,r1;
				pack2.unpack(r0,r1,helper_.leftRightSuper(ptr).left().
				             permutation(r));
				bool odd = superOddElectrons(t,ptr);
				odd ^= helper_.leftRightSuper(ptr).right().signs()[eta];
				const RealType sign = (odd && fermionSign ==
				                       ProgramGlobals::FermionOrBosonEnum::FERMION) ? -1.0
				                                                                    : 1.0;

				for (int k=Acrs.getRowPtr(r0);k<Acrs.getRowPtr(r0+1);k++) {
					SizeType r0prime = Acrs.getCol(k);
					for (int k2 = Bcrs.getRowPtr(eta);
					     k2<Bcrs.getRowPtr(eta+1);k2++) {
						SizeType eta2 = Bcrs.getCol(k2);
						SizeType rprime = helper_.leftRightSuper(ptr).left().
						        permutationInverse(r0prime+r1*ni);
						SizeType t2 = helper_.leftRightSuper(ptr).super().
						        permutationInverse(rprime+eta2*leftSize);
						if (t2<offset || t2>=total) continue;
						sum += Acrs.getValue(k)*Bcrs.getValue(k2)*
						        PsimagLite::conj(vec1.slowAccess(t))*
						        vec2.slowAccess(t2)*sign;
					}
				}
			}
		}

		return resultDivided(sum,vec1);
	}

	FieldType brLftCrnrEnviron_(const SparseMatrixType& Acrs,
	                            const SparseMatrixType& Bcrs,
	                            ProgramGlobals::FermionOrBosonEnum fOrB,
	                            const VectorWithOffsetType& vec1,
	                            const VectorWithOffsetType& vec2,
	                            SizeType ptr) const
	{
		const int fermionSign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;
		int signRight = fermionSignBasis(fermionSign,
		                                 helper_.leftRightSuper(ptr).right());
		FieldType sum = 0;
		SizeType ni = Bcrs.rows();
		SizeType leftSize = helper_.leftRightSuper(ptr).left().size();

		// some sanity checks:
		if (vec1.size() != vec2.size() ||
		        vec1.size()!=helper_.leftRightSuper(ptr).super().size())
			err("Observe::brLftCrnrEnviron_(...)\n");
		if (helper_.leftRightSuper(ptr).right().size()/Bcrs.rows() != Acrs.rows())
			err("Observe::brLftCrnrEnviron_(...)\n");

		// ok, we're ready for the main course:
		PackIndicesType pack1(helper_.leftRightSuper(ptr).left().size());
		PackIndicesType pack2(ni);

		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack1.unpack(eta,r,helper_.leftRightSuper(ptr).super().
				             permutation(t));
				SizeType r0,r1;
				pack2.unpack(r0,r1,helper_.leftRightSuper(ptr).right().permutation(r));
				const RealType sign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON)
				        ? 1 : helper_.signsOneSite(r0) * signRight;

				for (int k=Acrs.getRowPtr(r1);k<Acrs.getRowPtr(r1+1);k++) {
					SizeType r1prime = Acrs.getCol(k);
					for (int k2 = Bcrs.getRowPtr(eta);
					     k2<Bcrs.getRowPtr(eta+1);k2++) {
						SizeType eta2 = Bcrs.getCol(k2);
						SizeType rprime = helper_.leftRightSuper(ptr).right().
						        permutationInverse(r0+r1prime*ni);
						SizeType t2 = helper_.leftRightSuper(ptr).super().
						        permutationInverse(eta2+rprime*leftSize);
						if (t2<offset || t2>=total) continue;
						sum += PsimagLite::conj(Acrs.getValue(k))*Bcrs.getValue(k2)*
						        PsimagLite::conj(vec1.slowAccess(t))*
						        vec2.slowAccess(t2)*sign;
					}
				}
			}
		}

		return resultDivided(sum,vec1);
	}

	FieldType bracketRightCorner_(const SparseMatrixType& A1,
	                              const SparseMatrixType& A2,
	                              const SparseMatrixType& B,
	                              ProgramGlobals::FermionOrBosonEnum fOrB,
	                              const VectorWithOffsetType& vec1,
	                              const VectorWithOffsetType& vec2,
	                              SizeType ptr) const
	{
		if (helper_.direction(ptr) != ProgramGlobals::DirectionEnum::EXPAND_SYSTEM)
			return 0;

		const int fermionSign = (fOrB == ProgramGlobals::FermionOrBosonEnum::BOSON) ? 1 : -1;

		SparseMatrixType A1crs(A1);
		SparseMatrixType A2crs(A2);
		SparseMatrixType Bcrs(B);
		FieldType sum=0;
		SizeType ni = helper_.leftRightSuper(ptr).left().size()/Bcrs.rows();
		SizeType leftSize = helper_.leftRightSuper(ptr).left().size();

		// some sanity checks:
		assert(vec1.size()==vec2.size());

		if (vec1.size()==0) return 0;

		assert(vec1.size()==helper_.leftRightSuper(ptr).super().size());
		assert(ni==A1crs.rows());
		assert(Bcrs.rows()==A2crs.rows());

		// ok, we're ready for the main course:
		PackIndicesType pack1(helper_.leftRightSuper(ptr).left().size());
		PackIndicesType pack2(ni);

		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack1.unpack(r,
				             eta,
				             helper_.leftRightSuper(ptr).super().permutation(t));
				SizeType r0,r1;
				pack2.unpack(r0,
				             r1,
				             helper_.leftRightSuper(ptr).left().permutation(r));
				RealType sign = helper_.leftRightSuper(ptr).right().
				        fermionicSign(r1,fermionSign);

				for (int k1=A1crs.getRowPtr(r0);k1<A1crs.getRowPtr(r0+1);k1++) {
					SizeType r0prime = A1crs.getCol(k1);
					for (int k2=A2crs.getRowPtr(r1);k2<A2crs.getRowPtr(r1+1);k2++) {
						SizeType r1prime = A2crs.getCol(k2);
						for (int k3 = Bcrs.getRowPtr(eta);k3<Bcrs.getRowPtr(eta+1);k3++) {
							SizeType eta2 = Bcrs.getCol(k3);
							SizeType rprime = helper_.leftRightSuper(ptr).left().
							        permutationInverse(r0prime+r1prime*ni);
							SizeType t2 = helper_.leftRightSuper(ptr).super().
							        permutationInverse(rprime+eta2*leftSize);
							if (t2<offset || t2>=total) continue;
							sum +=  A1crs.getValue(k1)*A2crs.getValue(k2)*
							        Bcrs.getValue(k3)*
							        PsimagLite::conj(vec1.slowAccess(t))*
							        vec2.slowAccess(t2)*sign;
						}
					}
				}
			}
		}

		return resultDivided(sum,vec1);
	}

	FieldType resultDivided(FieldType sum, const VectorWithOffsetType& vec) const
	{
		FieldType tmp = vec*vec;
		RealType norma2 = PsimagLite::real(tmp);
		assert(fabs(norma2)>1e-10 && fabs(PsimagLite::imag(tmp))<1e-6);
		return sum/norma2;
	}

	const ObserverHelperType& helper_;
};  //class CorrelationsSkeleton
} // namespace Dmrg

/*@}*/
#endif

