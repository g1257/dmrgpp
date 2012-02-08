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
	typedef std::vector<ComplexOrRealType> VectorType;
	typedef SparseVector<typename VectorType::value_type> SparseVectorType;
	typedef ReflectionBasis<RealType,SparseMatrixType> ReflectionBasisType;

public:

	void update(const SparseMatrixType& sSector)
	{
		ReflectionBasisType reflectionBasis(sSector);

		computeTransform(Q1_,reflectionBasis,1.0);
		computeTransform(Qm_,reflectionBasis,-1.0);

	}

	void transform(SparseMatrixType& dest1,
		       SparseMatrixType& dest2,
		       const SparseMatrixType& src) const
	{

	}

	void setGs(VectorType& gs,const VectorType& v,size_t rank,size_t offset) const
	{
		VectorType gstmp(rank,0);

		for (size_t i=0;i<v.size();i++) {
			gstmp[i+offset]=v[i];
		}
		SparseMatrixType rT;
		//transposeConjugate(rT,transform_);
		multiply(gs,rT,gstmp);
	}

private:

	void computeTransform(SparseMatrixType& Q1,
			      const ReflectionBasisType& reflectionBasis,
			      const RealType& sector)
	{
		const SparseMatrixType& R1 = reflectionBasis.R(sector);
		SparseMatrixType R1Inverse;
		reflectionBasis.inverseTriangular(R1Inverse,R1);
		SparseMatrixType T1;

		buildT1(T1,reflectionBasis,sector);
		multiply(Q1,T1,R1Inverse);
	}

	void buildT1(SparseMatrixType& T1,
		     const ReflectionBasisType& reflectionBasis,
		     const RealType& sector) const
	{
		const std::vector<size_t>& ipPosOrNeg = reflectionBasis.ipPosOrNeg(sector);
		const SparseMatrixType& reflection = reflectionBasis.reflection();
		T1.resize(reflection.rank());
		size_t counter = 0;
		for (size_t i=0;i<reflection.rank();i++) {
			T1.setRow(i,counter);
			for (int k = reflection.getRowPtr(i);k<reflection.getRowPtr(i+1);k++) {
				size_t col = reflection.getCol(k);
				ComplexOrRealType val = reflection.getValue(k);
				if (col==i) {
					val += sector;
				}
				val *= sector;
				T1.pushCol(ipPosOrNeg[col]);
				T1.pushValue(val);
				counter++;
			}
		}
		T1.setRow(T1.rank(),counter);
		T1.checkValidity();
	}

	//	void checkTransform(const SparseMatrixType& sSector)
	//	{
	//		SparseMatrixType rT;
	//		transposeConjugate(rT,transform_);
	//		SparseMatrixType tmp3;
	//		multiply(tmp3,transform_,rT);

	////		printFullMatrix(rT,"Transform");

	//		SparseMatrixType tmp4;
	//		multiply(tmp3,sSector,rT);
	//		multiply(tmp4,transform_,tmp3);
	////		printFullMatrix(tmp4,"R S R^\\dagger");
	//	}

//	void split(SparseMatrixType& matrixA,SparseMatrixType& matrixB,const SparseMatrixType& matrix) const
//	{
//		size_t counter = 0;
//		matrixA.resize(plusSector_);
//		for (size_t i=0;i<plusSector_;i++) {
//			matrixA.setRow(i,counter);
//			for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
//				size_t col = matrix.getCol(k);
//				ComplexOrRealType val = matrix.getValue(k);
//				if (col<plusSector_) {
//					matrixA.pushCol(col);
//					matrixA.pushValue(val);
//					counter++;
//					continue;
//				}
//				if (!isAlmostZero(val,1e-5)) {
//					std::string s(__FILE__);
//					s += " Hamiltonian has no reflection symmetry.";
//					throw std::runtime_error(s.c_str());
//				}
//			}
//		}
//		matrixA.setRow(plusSector_,counter);

//		size_t rank = matrix.rank();
//		size_t minusSector=rank-plusSector_;
//		matrixB.resize(minusSector);
//		counter=0;
//		for (size_t i=plusSector_;i<rank;i++) {
//			matrixB.setRow(i-plusSector_,counter);
//			for (int k=matrix.getRowPtr(i);k<matrix.getRowPtr(i+1);k++) {
//				size_t col = matrix.getCol(k);
//				ComplexOrRealType val = matrix.getValue(k);
//				if (col>=plusSector_) {
//					matrixB.pushCol(col-plusSector_);
//					matrixB.pushValue(val);
//					counter++;
//					continue;
//				}
//				if (!isAlmostZero(val,1e-5)) {
//					std::string s(__FILE__);
//					s += " Hamiltonian has no reflection symmetry.";
//					throw std::runtime_error(s.c_str());
//				}
//			}
//		}
//		matrixB.setRow(minusSector,counter);
//	}

	SparseMatrixType Q1_,Qm_;

}; // class ReflectionTransform

} // namespace Dmrg 

/*@}*/
#endif // reflectionTRANSFORM_H
