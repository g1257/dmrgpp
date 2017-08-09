/*
Copyright (c)  2008-2013, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 3.0]
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
#include "Profiling.h"
#include "ApplyOperatorLocal.h"

namespace Dmrg {

// Don't add functions to this class
template<typename FieldType>
struct CorrelationData {
	SizeType ni;
	SizeType nj;
	typename PsimagLite::Vector<PsimagLite::Matrix<FieldType> >::Type correlationVector;
	PsimagLite::SparseVector<FieldType> wavefunction;
};

// Companion function:
template<typename FieldType>
std::ostream& operator<<(std::ostream& os,CorrelationData<FieldType>& c)
{
	os<<"ni="<<c.ni<<"\n";
	os<<"nj="<<c.nj<<"\n";
	os<<c.correlationVector;
	os<<c.wavefunction;
	return os;
}

template<typename ObserverHelperType_,typename ModelType>
class CorrelationsSkeleton {
	typedef SizeType IndexType;
	typedef PsimagLite::PackIndices PackIndicesType;

public:
	typedef ObserverHelperType_ ObserverHelperType;
	typedef typename ObserverHelperType::MatrixType MatrixType;
	typedef typename ObserverHelperType::VectorType VectorType ;
	typedef typename ObserverHelperType::VectorWithOffsetType VectorWithOffsetType;
	typedef typename ObserverHelperType::BasisWithOperatorsType BasisWithOperatorsType ;
	typedef typename BasisWithOperatorsType::RealType RealType;
	typedef typename BasisWithOperatorsType::BasisType BasisType;
	typedef typename ObserverHelperType::FermionSignType FermionSignType;

	typedef typename VectorType::value_type FieldType;
	typedef PsimagLite::Profiling ProfilingType;
	typedef typename BasisWithOperatorsType::OperatorType OperatorType;
	typedef PsimagLite::CrsMatrix<FieldType> SparseMatrixType;

	enum {GROW_RIGHT,GROW_LEFT};

	enum {DIAGONAL,NON_DIAGONAL};

	enum {LEFT_BRAKET=ObserverHelperType::LEFT_BRAKET,
		  RIGHT_BRAKET=ObserverHelperType::RIGHT_BRAKET};
	static const SizeType EXPAND_SYSTEM = ProgramGlobals::EXPAND_SYSTEM;
	static const SizeType EXPAND_ENVIRON = ProgramGlobals::EXPAND_ENVIRON;

	CorrelationsSkeleton(ObserverHelperType& helper,
	                     const ModelType&,
	                     bool verbose = false)
	    : helper_(helper),verbose_(verbose)
	{}

	SizeType numberOfSites(SizeType threadId) const
	{
		return helper_.leftRightSuper(threadId).sites();
	}

	//! i can be zero here!!
	void growDirectly(SparseMatrixType& Odest,
	                  const SparseMatrixType& Osrc,
	                  SizeType i,
	                  int fermionicSign,
	                  SizeType ns,
	                  bool transform,
	                  SizeType threadId)
	{
		Odest =Osrc;
		// from 0 --> i
		int nt=i-1;
		if (nt<0) nt=0;

		for (SizeType s=nt;s<ns;s++) {
			helper_.setPointer(threadId,s);
			SizeType growOption = growthDirection(s,nt,i,threadId);
			SparseMatrixType Onew(helper_.columns(threadId),helper_.columns(threadId));

			fluffUp(Onew,Odest,fermionicSign,growOption,false,threadId);
			if (!transform && s == ns-1) {
				Odest = Onew;
				continue;
			}

			helper_.transform(Odest,Onew,threadId);
		}
	}

	SizeType growthDirection(SizeType s,int nt,SizeType i,SizeType threadId) const
	{
		SizeType dir = helper_.direction(threadId);
		SizeType growOption = (dir==EXPAND_SYSTEM) ? GROW_RIGHT : GROW_LEFT;

		if (s==SizeType(nt)) {
			growOption = (dir==EXPAND_SYSTEM) ? GROW_LEFT : GROW_RIGHT;
			if (i==0) growOption = (dir==EXPAND_SYSTEM) ? GROW_RIGHT : GROW_LEFT;
		}
		return growOption;
	}

	// Perfomance critical:
	void fluffUp(SparseMatrixType& ret2,
	             const SparseMatrixType& O,
	             int fermionicSign,
	             int growOption,
	             bool transform,
	             SizeType threadId)
	{
		if (helper_.direction(threadId)==EXPAND_SYSTEM) {
			fluffUpSystem(ret2,O,fermionicSign,growOption,transform,threadId);
			return;
		}
		fluffUpEnviron(ret2,O,fermionicSign,growOption,transform,threadId);
	}

	void dmrgMultiply(SparseMatrixType& result,
	                  const SparseMatrixType& O1,
	                  const SparseMatrixType& O2,
	                  int fermionicSign,
	                  SizeType ns,
	                  SizeType threadId)
	{
		if (helper_.direction(threadId)==EXPAND_SYSTEM) {
			dmrgMultiplySystem(result,O1,O2,fermionicSign,ns,threadId);
			return;
		}
		dmrgMultiplyEnviron(result,O1,O2,fermionicSign,ns,threadId);
	}

	void createWithModification(SparseMatrixType& Om,const SparseMatrixType& O,char mod)
	{
		if (mod == 'n' || mod == 'N') {
			Om = O;
			return;
		}

		transposeConjugate(Om,O);
	}

	FieldType bracket(const SparseMatrixType& A,int fermionicSign,SizeType threadId)
	{
		try {
			const VectorWithOffsetType& src1 =
			        helper_.getVectorFromBracketId(LEFT_BRAKET,threadId);
			const VectorWithOffsetType& src2 =
			        helper_.getVectorFromBracketId(RIGHT_BRAKET,threadId);

			return bracket_(A,src1,src2,fermionicSign,threadId);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: CorrelationsSkeleton::bracket(...):";
			std::cerr<<" No data seen yet\n";
			return 0;
		}
	}

	FieldType bracketRightCorner(const SparseMatrixType& A,
	                             const SparseMatrixType& B,
	                             int fermionSign,
	                             SizeType threadId)
	{
		try {
			const VectorWithOffsetType& src1 =
			        helper_.getVectorFromBracketId(LEFT_BRAKET,threadId);
			const VectorWithOffsetType& src2 =
			        helper_.getVectorFromBracketId(RIGHT_BRAKET,threadId);
			return bracketRightCorner_(A,B,fermionSign,src1,src2,threadId);
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
	                             int fermionSign,
	                             SizeType threadId)
	{
		try {
			const VectorWithOffsetType& src1 =
			        helper_.getVectorFromBracketId(LEFT_BRAKET,threadId);
			const VectorWithOffsetType& src2 =
			        helper_.getVectorFromBracketId(RIGHT_BRAKET,threadId);
			return bracketRightCorner_(A,B,C,fermionSign,src1,src2,threadId);
		} catch (std::exception& e) {
			std::cerr<<"CAUGHT: "<<e.what();
			std::cerr<<"WARNING: CorrelationsSkeleton::bracketRightCornerABC(...):";
			std::cerr<<" No data seen yet\n";
			return 0;
		}
	}

private:

	void dmrgMultiplySystem(SparseMatrixType& result,
	                        const SparseMatrixType& O1,
	                        const SparseMatrixType& O2,
	                        int fermionicSign,
	                        SizeType ns,
	                        SizeType threadId)
	{
		SizeType ni=O1.rows();

		helper_.setPointer(threadId,ns);
		SizeType sprime = helper_.leftRightSuper(threadId).left().size(); //ni*nj;
		result.resize(sprime,sprime);

		if (helper_.leftRightSuper(threadId).left().size()!=sprime) {
			std::cerr<<"WARNING: "<<helper_.leftRightSuper(threadId).left().size();
			std::cerr<<"!="<<sprime<<"\n";
			throw PsimagLite::RuntimeError("problem in dmrgMultiply\n");
		}

		PsimagLite::Vector<SizeType>::Type col(sprime,0);
		typename PsimagLite::Vector<FieldType>::Type value(sprime,0);

		PackIndicesType pack(ni);

		SizeType counter = 0;
		for (SizeType r=0;r<sprime;r++) {
			SizeType e,u;
			pack.unpack(e,u,helper_.leftRightSuper(threadId).left().permutation(r));
			RealType f = helper_.fermionicSignLeft(threadId)(e,fermionicSign);
			result.setRow(r,counter);
			for (int k=O1.getRowPtr(e);k<O1.getRowPtr(e+1);k++) {
				SizeType e2 = O1.getCol(k);
				for (int k2=O2.getRowPtr(u);k2<O2.getRowPtr(u+1);k2++) {
					SizeType u2 = O2.getCol(k2);
					SizeType r2 = helper_.leftRightSuper(threadId).left().
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
	                         const SparseMatrixType& O1,
	                         const SparseMatrixType& O2,
	                         int fermionicSign,
	                         SizeType ns,
	                         SizeType threadId)
	{
		SizeType nj=O2.rows();

		helper_.setPointer(threadId,ns);
		SizeType eprime = helper_.leftRightSuper(threadId).right().size(); //ni*nj;
		result.resize(eprime,eprime);

		if (helper_.leftRightSuper(threadId).right().size()!=eprime) {
			std::cerr<<"WARNING: "<<helper_.leftRightSuper(threadId).right().size();
			std::cerr<<"!="<<eprime<<"\n";
			throw PsimagLite::RuntimeError("problem in dmrgMultiply\n");
		}

		PsimagLite::Vector<SizeType>::Type col(eprime,0);
		typename PsimagLite::Vector<FieldType>::Type value(eprime,0);

		PackIndicesType pack(nj);

		SizeType counter = 0;
		for (SizeType r=0;r<eprime;r++) {
			result.setRow(r,counter);
			SizeType e,u;

			pack.unpack(e,u,helper_.leftRightSuper(threadId).right().permutation(r));
			SizeType nx0 = helper_.leftRightSuper(threadId).right().
			        electrons(BasisType::AFTER_TRANSFORM);
			RealType f = (nx0 & 1) ? fermionicSign : 1;

			for (int k=O2.getRowPtr(e);k<O2.getRowPtr(e+1);k++) {
				SizeType e2 = O2.getCol(k);
				for (int k2=O1.getRowPtr(u);k2<O1.getRowPtr(u+1);k2++) {
					SizeType u2 = O1.getCol(k2);
					SizeType r2 = helper_.leftRightSuper(threadId).right().
					        permutationInverse(e2 + u2*nj);
					assert(r2<eprime);
					col[r2] = 1;
					value[r2] += O2.getValue(k)*O1.getValue(k2)*f;
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

	// Perfomance critical:
	void fluffUpSystem(SparseMatrixType& ret2,
	                   const SparseMatrixType& O,
	                   int fermionicSign,
	                   int growOption,
	                   bool transform,
	                   SizeType threadId)
	{
		SizeType n =helper_.leftRightSuper(threadId).left().size();

		MatrixType ret(n,n);

		for (SizeType e=0;e<n;e++) {
			for (SizeType e2=0;e2<n;e2++) {
				ret(e,e2) = fluffUpSystem_(
				            O,e,e2,fermionicSign,growOption,threadId);
			}
		}

		if (transform) {
			SparseMatrixType ret3(ret);
			helper_.transform(ret2,ret3,threadId);
			return;
		}

		fullMatrixToCrsMatrix(ret2,ret);
	}

	// Perfomance critical:
	void fluffUpEnviron(SparseMatrixType& ret2,
	                    const SparseMatrixType& O,
	                    int fermionicSign,
	                    int growOption,
	                    bool transform,
	                    SizeType threadId)
	{
		SizeType n =helper_.leftRightSuper(threadId).right().size();

		MatrixType ret(n,n);
		for (SizeType e=0;e<n;e++) {
			for (SizeType e2=0;e2<n;e2++) {
				ret(e,e2) = fluffUpEnviron_(
				            O,e,e2,fermionicSign,growOption,threadId);
			}
		}
		if (transform) {
			SparseMatrixType ret3(ret);
			helper_.transform(ret2,ret3,threadId);
			return;
		}

		fullMatrixToCrsMatrix(ret2,ret);
	}

	// Perfomance critical:
	FieldType fluffUpSystem_(const SparseMatrixType& O,
	                         SizeType e,SizeType e2,
	                         int fermionicSign,
	                         int growOption,
	                         SizeType threadId)
	{
		SizeType n = O.rows();
		SizeType m = SizeType(helper_.leftRightSuper(threadId).left().size()/n);
		RealType sign = static_cast<RealType>(1.0);

		// Sperm[e] = i +k*n or e= k + i*m
		// Sperm[e2] = j+k*n or e2=k+j*m
		SizeType i,j,k,k2;
		if (growOption==GROW_RIGHT) {
			if (SizeType(helper_.leftRightSuper(threadId).left().permutation(e)/n)!=
			    SizeType(helper_.leftRightSuper(threadId).left().permutation(e2)/n))
				return 0;
			PackIndicesType pack(n);
			pack.unpack(i,k,helper_.leftRightSuper(threadId).left().permutation(e));
			pack.unpack(j,k2,helper_.leftRightSuper(threadId).left().permutation(e2));
		} else {
			if (SizeType(helper_.leftRightSuper(threadId).left().permutation(e)%m)!=
			    SizeType(helper_.leftRightSuper(threadId).left().permutation(e2)%m))
				return 0;
			PackIndicesType pack(m);
			pack.unpack(k,i,helper_.leftRightSuper(threadId).left().permutation(e));
			pack.unpack(k2,j,helper_.leftRightSuper(threadId).left().permutation(e2));
			sign = helper_.fermionicSignLeft(threadId)(k,fermionicSign);
		}
		if (k!=k2) return 0;
		return O.element(i,j)*sign;
	}

	// Perfomance critical:
	FieldType fluffUpEnviron_(const SparseMatrixType& O,
	                          SizeType e,SizeType e2,
	                          int fermionicSign,
	                          int growOption,
	                          SizeType threadId)
	{
		SizeType n = O.rows();
		SizeType m = SizeType(helper_.leftRightSuper(threadId).right().size()/n);
		RealType sign = 1;

		// Eperm[e] = i +k*n or e= k + i*m
		// Eperm[e2] = j+k*n or e2=k+j*m

		SizeType i,j,k,k2;

		if (growOption==GROW_RIGHT) {
			PackIndicesType pack(n);
			pack.unpack(i,k,helper_.leftRightSuper(threadId).right().permutation(e));
			pack.unpack(j,k2,helper_.leftRightSuper(threadId).right().permutation(e2));
			SizeType nx0 = helper_.leftRightSuper(threadId).left().
			        electrons(BasisType::AFTER_TRANSFORM);
			sign = (nx0 & 1) ? fermionicSign : 1;
		} else {
			PackIndicesType pack(m);
			pack.unpack(k,i,helper_.leftRightSuper(threadId).right().permutation(e));
			pack.unpack(k2,j,helper_.leftRightSuper(threadId).right().permutation(e2));
			SizeType nx0 = helper_.leftRightSuper(threadId).super().
			        electrons(BasisType::AFTER_TRANSFORM);
			sign = (nx0 & 1) ?  fermionicSign : 1;
		}
		if (k!=k2) return 0;
		return O.element(i,j)*sign;
	}

	FieldType bracket_(const SparseMatrixType& A,
	                   const VectorWithOffsetType& vec1,
	                   const VectorWithOffsetType& vec2,
	                   int fermionicSign,
	                   SizeType threadId)
	{
		if (verbose_)
			std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

		if (vec1.size()!=helper_.leftRightSuper(threadId).super().size() ||
		    vec1.size()!=vec2.size())
			throw PsimagLite::RuntimeError("CorrelationsSkeleton::bracket_(...): Error\n");

		if (helper_.direction(threadId)==EXPAND_SYSTEM) {
			return bracketSystem_(A,vec1,vec2,threadId);
		}

		return bracketEnviron_(A,vec1,vec2,fermionicSign,threadId);
	}

	FieldType bracketSystem_(const SparseMatrixType& A,
	                        const VectorWithOffsetType& vec1,
	                        const VectorWithOffsetType& vec2,
	                        SizeType threadId)
	{
		FieldType sum=0;
		PackIndicesType pack(helper_.leftRightSuper(threadId).left().size());
		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack.unpack(r,eta,helper_.leftRightSuper(threadId).super().
				            permutation(t));
				for (int k=A.getRowPtr(r);k<A.getRowPtr(r+1);k++) {
					SizeType r2 = A.getCol(k);
					SizeType t2 = helper_.leftRightSuper(threadId).super().
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
	                         int fermionicSign,
	                         SizeType threadId)
	{
		FieldType sum=0;
		PackIndicesType pack(helper_.leftRightSuper(threadId).left().size());
		SizeType leftSize = helper_.leftRightSuper(threadId).left().size();

		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack.unpack(r,eta,helper_.leftRightSuper(threadId).super().
				            permutation(t));
				if (eta>=A.rows()) throw PsimagLite::RuntimeError("Error\n");
				SizeType nx0 = helper_.leftRightSuper(threadId).left().
				        electrons(BasisType::AFTER_TRANSFORM);
				RealType sign = (nx0 & 1) ? fermionicSign : 1;
				for (int k=A.getRowPtr(eta);k<A.getRowPtr(eta+1);k++) {
					SizeType eta2 = A.getCol(k);
					SizeType t2 = helper_.leftRightSuper(threadId).super().
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
	                             int fermionSign,
	                             const VectorWithOffsetType& vec1,
	                             const VectorWithOffsetType& vec2,
	                             SizeType threadId)
	{
		if (helper_.direction(threadId)==EXPAND_SYSTEM)
			return brRghtCrnrSystem_(A,B,fermionSign,vec1,vec2,threadId);
		return brLftCrnrEnviron_(A,B,fermionSign,vec1,vec2,threadId);
	}

	SizeType superElectrons(SizeType t, SizeType threadId) const
	{
#if 0
		return helper_.leftRightSuper(threadId).super().electrons(t);
#else
		SizeType tmp = helper_.leftRightSuper(threadId).super().permutation(t);
		div_t mydiv = div(tmp,helper_.leftRightSuper(threadId).left().size());
		return helper_.leftRightSuper(threadId).right().electrons(mydiv.quot) +
		        helper_.leftRightSuper(threadId).left().electrons(mydiv.rem);
#endif

	}

	FieldType brRghtCrnrSystem_(const SparseMatrixType& Acrs,
	                           const SparseMatrixType& Bcrs,
	                           int fermionSign,
	                           const VectorWithOffsetType& vec1,
	                           const VectorWithOffsetType& vec2,
	                           SizeType threadId)
	{
		if (verbose_)
			std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

		FieldType sum=0;
		SizeType leftSize = helper_.leftRightSuper(threadId).left().size();
		SizeType ni = helper_.leftRightSuper(threadId).left().size()/Bcrs.rows();

		// some sanity checks:
		if (vec1.size()!=vec2.size() || vec1.size()!=
		    helper_.leftRightSuper(threadId).super().size())
			throw PsimagLite::RuntimeError("Observe::brRghtCrnrSystem_(...)\n");
		if (ni!=Acrs.rows())
			throw PsimagLite::RuntimeError("Observe::brRghtCrnrSystem_(...)\n");

		// ok, we're ready for the main course:
		PackIndicesType pack1(helper_.leftRightSuper(threadId).left().size());
		PackIndicesType pack2(ni);
		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack1.unpack(r,eta,helper_.leftRightSuper(threadId).super().
				             permutation(t));
				SizeType r0,r1;
				pack2.unpack(r0,r1,helper_.leftRightSuper(threadId).left().
				             permutation(r));
				SizeType electrons = superElectrons(t,threadId);
				electrons -= helper_.leftRightSuper(threadId).right().electrons(eta);
				RealType sign = (electrons & 1) ? fermionSign : 1.0;

				for (int k=Acrs.getRowPtr(r0);k<Acrs.getRowPtr(r0+1);k++) {
					SizeType r0prime = Acrs.getCol(k);
					for (int k2 = Bcrs.getRowPtr(eta);
					     k2<Bcrs.getRowPtr(eta+1);k2++) {
						SizeType eta2 = Bcrs.getCol(k2);
						SizeType rprime = helper_.leftRightSuper(threadId).left().
						        permutationInverse(r0prime+r1*ni);
						SizeType t2 = helper_.leftRightSuper(threadId).super().
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
	                           int fermionSign,
	                           const VectorWithOffsetType& vec1,
	                           const VectorWithOffsetType& vec2,
	                           SizeType threadId)
	{
		if (verbose_)
			std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

		FieldType sum=0;
		SizeType ni = Bcrs.rows();
		SizeType leftSize = helper_.leftRightSuper(threadId).left().size();

		// some sanity checks:
		if (vec1.size()!=vec2.size() ||
		    vec1.size()!=helper_.leftRightSuper(threadId).super().size())
			throw PsimagLite::RuntimeError("Observe::brLftCrnrEnviron_(...)\n");
		if (helper_.leftRightSuper(threadId).right().size()/Bcrs.rows()!=Acrs.rows())
			throw PsimagLite::RuntimeError("Observe::bracketRightCorner_(...)\n");

		// ok, we're ready for the main course:
		PackIndicesType pack1(helper_.leftRightSuper(threadId).left().size());
		PackIndicesType pack2(ni);

		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack1.unpack(eta,r,helper_.leftRightSuper(threadId).super().
				             permutation(t));
				SizeType r0,r1;
				pack2.unpack(r0,r1,helper_.leftRightSuper(threadId).right().permutation(r));
				SizeType electrons = helper_.leftRightSuper(threadId).left().electrons(eta);
				RealType sign = (electrons & 1) ? fermionSign : 1.0;

				for (int k=Acrs.getRowPtr(r1);k<Acrs.getRowPtr(r1+1);k++) {
					SizeType r1prime = Acrs.getCol(k);
					for (int k2 = Bcrs.getRowPtr(eta);
					     k2<Bcrs.getRowPtr(eta+1);k2++) {
						SizeType eta2 = Bcrs.getCol(k2);
						SizeType rprime = helper_.leftRightSuper(threadId).right().
						        permutationInverse(r0+r1prime*ni);
						SizeType t2 = helper_.leftRightSuper(threadId).super().
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
	                             int fermionSign,
	                             const VectorWithOffsetType& vec1,
	                             const VectorWithOffsetType& vec2,
	                             SizeType threadId)
	{
		if (helper_.direction(threadId)!=EXPAND_SYSTEM) return 0;

		if (verbose_)
			std::cerr<<"SE.size="<<helper_.leftRightSuper(threadId).super().size()<<"\n";

		SparseMatrixType A1crs(A1);
		SparseMatrixType A2crs(A2);
		SparseMatrixType Bcrs(B);
		FieldType sum=0;
		SizeType ni = helper_.leftRightSuper(threadId).left().size()/Bcrs.rows();
		SizeType leftSize = helper_.leftRightSuper(threadId).left().size();

		// some sanity checks:
		assert(vec1.size()==vec2.size());

		if (vec1.size()==0) return 0;

		assert(vec1.size()==helper_.leftRightSuper(threadId).super().size());
		assert(ni==A1crs.rows());
		assert(Bcrs.rows()==A2crs.rows());

		// ok, we're ready for the main course:
		PackIndicesType pack1(helper_.leftRightSuper(threadId).left().size());
		PackIndicesType pack2(ni);

		for (SizeType x=0;x<vec1.sectors();x++) {
			SizeType sector = vec1.sector(x);
			SizeType offset = vec1.offset(sector);
			SizeType total = offset + vec1.effectiveSize(sector);
			for (SizeType t=offset;t<total;t++) {
				SizeType eta,r;

				pack1.unpack(r,
				             eta,
				             helper_.leftRightSuper(threadId).super().permutation(t));
				SizeType r0,r1;
				pack2.unpack(r0,
				             r1,
				             helper_.leftRightSuper(threadId).left().permutation(r));
				RealType sign = helper_.leftRightSuper(threadId).right().
				        fermionicSign(r1,fermionSign);

				for (int k1=A1crs.getRowPtr(r0);k1<A1crs.getRowPtr(r0+1);k1++) {
					SizeType r0prime = A1crs.getCol(k1);
					for (int k2=A2crs.getRowPtr(r1);k2<A2crs.getRowPtr(r1+1);k2++) {
						SizeType r1prime = A2crs.getCol(k2);
						for (int k3 = Bcrs.getRowPtr(eta);k3<Bcrs.getRowPtr(eta+1);k3++) {
							SizeType eta2 = Bcrs.getCol(k3);
							SizeType rprime = helper_.leftRightSuper(threadId).left().
							        permutationInverse(r0prime+r1prime*ni);
							SizeType t2 = helper_.leftRightSuper(threadId).super().
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

	ObserverHelperType& helper_; //<-- NB: We are not the owner
	bool verbose_;
};  //class CorrelationsSkeleton
} // namespace Dmrg

/*@}*/
#endif

