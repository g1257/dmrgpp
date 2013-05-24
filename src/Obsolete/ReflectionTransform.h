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

/*! \file ReflectionTransform
 *
 *
 */
#ifndef reflectionTRANSFORM_H
#define reflectionTRANSFORM_H

#include "PackIndices.h" // in PsimagLite
#include "Matrix.h"
#include "ProgressIndicator.h"
#include "LAPACK.h"
#include "Sort.h"
#include "ReflectionBasis.h"

namespace Dmrg {

template<typename RealType,typename SparseMatrixType>
class ReflectionTransform {

	typedef typename SparseMatrixType::value_type ComplexOrRealType;
	typedef typename PsimagLite::Vector<ComplexOrRealType>::Type VectorType;
	typedef SparseVector<typename VectorType::value_type> SparseVectorType;
	typedef ReflectionBasis<RealType,SparseMatrixType> ReflectionBasisType;

public:

	ReflectionTransform(bool idebug)
	: idebug_(idebug)
	{}

	void update(const SparseMatrixType& sSector)
	{
		ReflectionBasisType reflectionBasis(sSector,idebug_);
		plusSector_ = reflectionBasis.R(1.0).rank();
		computeTransform(Q1_,reflectionBasis,1.0);
		computeTransform(Qm_,reflectionBasis,-1.0);
//		SparseMatrixType Q;
		//computeFullQ(Q,Q1_,Qm_);
//		split(Q);
		if (!idebug_) return;
		printFullMatrix(Q1_,"Q1");
		printFullMatrix(Qm_,"Qm");
	}

	void transform(SparseMatrixType& dest1,
		       SparseMatrixType& destm,
		       const SparseMatrixType& H) const
	{
		SparseMatrixType HQ1,HQm;
		multiply(HQ1,H,Q1_);

		multiply(HQm,H,Qm_);
		if (idebug_) {
			printFullMatrix(H,"OriginalHamiltonian");
			printFullMatrix(HQm,"HQm");
			printFullMatrix(HQ1,"HQ1");
		}

		SparseMatrixType Q1t,Qmt;
		transposeConjugate(Q1t,Q1_);
		transposeConjugate(Qmt,Qm_);

		RealType norm1 = getNorm(Q1t,Q1_);
		RealType normm = getNorm(Qmt,Qm_);


		multiply(dest1,Q1t,HQ1);
		reshape(dest1,plusSector_);
		dest1 *= (1.0/norm1);
		assert(isHermitian(dest1));

		SizeType minusSector = H.rank()-plusSector_;

		multiply(destm,Qmt,HQm);
		reshape(destm,minusSector);
		destm *= (1.0/normm);
		assert(isHermitian(destm));

		if (idebug_) {
			std::cerr<<"norm1="<<norm1<<" normm="<<normm<<"\n";
			std::cerr<<"plusSector="<<plusSector_<<" minusSector="<<minusSector<<"\n";
			printFullMatrix(dest1,"dest1");
			printFullMatrix(destm,"destm");
		}

#ifndef NDEBUG
//		checkTransform(Qmt,HQ1);
//		checkTransform(Q1t,HQm);
//		SparseMatrixType A;
//		multiply(A,Q1t,Q1_);
//		bool b = isThePartialIdentity(A,plusSector_,1e-5);
//		if (!b) {
//			printFullMatrix(A,"A");
//			assert(b);
//		}

//		multiply(A,Qmt,Qm_);
//		b = isThePartialIdentity(A,A.rank()-plusSector_,1e-5);
//		if (!b) {
//			printFullMatrix(A,"A");
//			assert(b);
//		}
#endif
	}

	bool isThePartialIdentity(const SparseMatrixType& A,SizeType partialSize,const RealType& eps = 1e-6) const
	{
		for (SizeType i=0;i<partialSize;i++) {
			for (int k = A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				SizeType col = A.getCol(k);
				if (col>=partialSize) continue;
				ComplexOrRealType val = A.getValue(k);
				if (i==col && !isAlmostZero(val-1.0,eps)) return false;
				if (i!=col && !isAlmostZero(val,eps)) return false;
			}
		}
		return true;
	}

	void setGs(VectorType& gs,const VectorType& v,const RealType& sector) const
	{
		const SparseMatrixType& Q = (sector>0) ? Q1_ : Qm_;
		multiply(gs,Q,v);
		RealType norma = PsimagLite::norm(gs);
		gs /= norma;
	}

	template<typename SomeVectorType>
	void setInitState(const SomeVectorType& initVector,
			  SomeVectorType& initVector1,
			  SomeVectorType& initVector2) const
	{
		SizeType minusSector = initVector.size()-plusSector_;
		initVector1.resize(plusSector_);
		initVector2.resize(minusSector);
		for (SizeType i=0;i<initVector.size();i++) {
			if (i<plusSector_) initVector1[i]=initVector[i];
			else  initVector2[i-plusSector_]=initVector[i];
		}
	}

private:

	RealType getNorm(const SparseMatrixType& A,const SparseMatrixType& B) const
	{
		SparseMatrixType C;
		multiply(C,A,B);
		for (SizeType i=0;i<C.rank();i++) {
			for (int k=C.getRowPtr(i);k<C.getRowPtr(i+1);k++) {
				SizeType col = C.getCol(k);
				if (col==i) {
					return std::real(C.getValue(k));
				}
			}
		}
		assert(false);
		return 0;
	}

	void reshape(SparseMatrixType& A,SizeType n2) const
	{
		SparseMatrixType B(n2,n2);
		SizeType counter = 0;
		for (SizeType i=0;i<n2;i++) {
			B.setRow(i,counter);
			for (int k = A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				SizeType col = A.getCol(k);
				ComplexOrRealType val = A.getValue(k);
				if (col>=n2) {
					assert(isAlmostZero(val,1e-5));
					continue;
				}
				B.pushCol(col);
				B.pushValue(val);
				counter++;
			}
		}
#ifndef NDEBUG
		for (SizeType i=n2;i<A.rank();i++) {
			for (int k = A.getRowPtr(i);k<A.getRowPtr(i+1);k++) {
				ComplexOrRealType val = A.getValue(k);
				assert(isAlmostZero(val,1e-5));
			}
		}
#endif
		B.setRow(n2,counter);
		B.checkValidity();
		A = B;
	}

	void checkTransform(const SparseMatrixType& A,const SparseMatrixType& B) const
	{
		SparseMatrixType C;
		multiply(C,A,B);
		bool b = isZero(C,1e-5);
		if (b) return;
		printFullMatrix(A,"MatrixA");
		printFullMatrix(B,"MatrixB");
		printFullMatrix(C,"MatrixC");
		assert(b);
	}

	void computeTransform(SparseMatrixType& Q1,
			      const ReflectionBasisType& reflectionBasis,
			      const RealType& sector)
	{
		const SparseMatrixType& R1 = reflectionBasis.R(sector);
		SparseMatrixType R1Inverse;
		reflectionBasis.inverseTriangular(R1Inverse,R1,sector);

		SparseMatrixType T1;

		buildT1(T1,reflectionBasis,sector);
		bool strict = false; // matrices below have different ranks!!
		multiply(Q1,T1,R1Inverse,strict);
		if (!idebug_) return;
		printFullMatrix(R1Inverse,"R1Inverse");
		printFullMatrix(T1,"T1");
	}

	void computeFullQ(SparseMatrixType& Q,
			  const SparseMatrixType& Q1,
			  const SparseMatrixType& Qm) const
	{
		SizeType n = Q1.rank();
		typename PsimagLite::Vector<ComplexOrRealType>::Type sum(n,0.0);
		SizeType counter = 0;
		Q.resize(n);
		SizeType minusSector = n - plusSector_;
		for (SizeType i=0;i<n;i++) {
			Q.setRow(i,counter);
			// add Q1
			for (int k = Q1.getRowPtr(i);k<Q1.getRowPtr(i+1);k++) {
				SizeType col = Q1.getCol(k);
				if (col>=plusSector_) continue;
				ComplexOrRealType val =  Q1.getValue(k);
				Q.pushValue(val);
				Q.pushCol(col);
				sum[i] += std::conj(val)*val;
				counter++;
			}
			// add Qm
			for (int k = Qm.getRowPtr(i);k<Qm.getRowPtr(i+1);k++) {
				SizeType col = Qm.getCol(k);
				if (col>=minusSector) continue;
				ComplexOrRealType val =  Qm.getValue(k);
				Q.pushValue(val);
				Q.pushCol(Qm.getCol(k)+plusSector_);
				sum[i] += std::conj(val)*val;
				counter++;
			}
		}
		Q.setRow(Q.rank(),counter);
		// normalize
//		for (SizeType i=0;i<n;i++) {
//			if (isAlmostZero(sum[i],1e-10)) continue;
//			sum[i] = 1.0/sqrt(sum[i]);
//			for (int k = Q.getRowPtr(i);k<Q.getRowPtr(i+1);k++) {
//				Q.setValues(k,Q.getValue(k)*sum[i]);
//			}
//		}
		Q.checkValidity();
#ifndef NDEBUG
		SparseMatrixType Qt;
		transposeConjugate(Qt,Q);
		SparseMatrixType A;
		multiply(A,Qt,Q);
		if (!isTheIdentity(A,1e-5)) {
			printFullMatrix(Q,"Q");
			printFullMatrix(A,"A");
			assert(false);
		}
#endif
		if (!idebug_) return;
		printFullMatrix(Q,"Q");
	}

	void split(const SparseMatrixType& Q)
	{
		SizeType n = Q.rank();
		Q1_.resize(n);
		SizeType counter = 0;
		for (SizeType i=0;i<n;i++) {
			Q1_.setRow(i,counter);
			for (int k = Q.getRowPtr(i);k<Q.getRowPtr(i+1);k++) {
				SizeType col = Q.getCol(k);
				if (col>=plusSector_) continue;
				Q1_.pushCol(col);
				Q1_.pushValue(Q.getValue(k));
				counter++;
			}
		}
		Q1_.setRow(Q1_.rank(),counter);

		counter = 0;
		Qm_.resize(n);
		for (SizeType i=0;i<n;i++) {
			Qm_.setRow(i,counter);
			for (int k = Q.getRowPtr(i);k<Q.getRowPtr(i+1);k++) {
				SizeType col = Q.getCol(k);
				if (col<plusSector_) continue;
				Qm_.pushCol(col-plusSector_);
				Qm_.pushValue(Q.getValue(k));
				counter++;
			}
		}
		Qm_.setRow(Qm_.rank(),counter);
	}

	void buildT1(SparseMatrixType& T1final,
		     const ReflectionBasisType& reflectionBasis,
		     const RealType& sector) const
	{
		const typename PsimagLite::Vector<SizeType>::Type& ipPosOrNeg = reflectionBasis.ipPosOrNeg(sector);
		const SparseMatrixType& reflection = reflectionBasis.reflection();
		SizeType n = reflection.rank();

		SparseMatrixType T1(n,n);
		SizeType counter = 0;
		for (SizeType i=0;i<n;i++) {
			T1.setRow(i,counter);
			bool hasDiagonal = false;
			for (int k = reflection.getRowPtr(i);k<reflection.getRowPtr(i+1);k++) {
				SizeType col = reflection.getCol(k);
				ComplexOrRealType val = reflection.getValue(k);
				if (col==i) {
					val += sector;
					hasDiagonal=true;
				}
				val *= sector;
				if (isAlmostZero(val,1e-10)) continue;
				T1.pushCol(col);
				T1.pushValue(val);
				counter++;
			}
			if (!hasDiagonal) {
				T1.pushCol(i);
				T1.pushValue(1.0);
				counter++;
			}

		}
		T1.setRow(n,counter);
		T1.checkValidity();

		// permute columns now:
		typename PsimagLite::Vector<int>::Type inversePermutation(n,-1);
		for (SizeType i=0;i<ipPosOrNeg.size();i++)
			inversePermutation[ipPosOrNeg[i]]=i;

		T1final.resize(n);
		counter=0;
		typename PsimagLite::Vector<ComplexOrRealType>::Type sum(n,0.0);
		for (SizeType i=0;i<n;i++) {
			T1final.setRow(i,counter);
			for (int k = T1.getRowPtr(i);k<T1.getRowPtr(i+1);k++) {
				int col = inversePermutation[T1.getCol(k)];
				if (col<0) continue;
				ComplexOrRealType val = T1.getValue(k);
				if (isAlmostZero(val,1e-10)) continue;
				assert(SizeType(col)<ipPosOrNeg.size());
				T1final.pushCol(col);
				T1final.pushValue(val);
				sum[i] += std::conj(val)*val;
				counter++;
			}
		}
		T1final.setRow(n,counter);
		T1final.checkValidity();

		// normalize T1
//		for (SizeType i=0;i<n;i++) {
//			if (isAlmostZero(sum[i],1e-12)) continue;
//			sum[i] = 1.0/sqrt(sum[i]);

//			for (int k = T1final.getRowPtr(i);k<T1final.getRowPtr(i+1);k++)
//				T1final.setValues(k,T1final.getValue(k) * sum[i]);

//		}
	}

	bool idebug_;
	SizeType plusSector_;
	SparseMatrixType Q1_,Qm_;

}; // class ReflectionTransform

} // namespace Dmrg 

/*@}*/
#endif // reflectionTRANSFORM_H
