/*
Copyright (c) 2009-2017, UT-Battelle, LLC
All rights reserved

[DMRG++, Version 4.]
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

/*! \file DensityMatrixSu2.h
 *
 *
 *
 */
#ifndef DENSITY_MATRIX_SU2_H
#define DENSITY_MATRIX_SU2_H

#include "Matrix.h" // in PsimagLite
#include "AlmostEqual.h" // in PsimagLite
#include "BlockDiagonalMatrix.h"
#include "DensityMatrixBase.h"
#include "ProgramGlobals.h"
#include "DiagBlockDiagMatrix.h"

namespace Dmrg {
template<typename TargetingType>
class DensityMatrixSu2 : public DensityMatrixBase<TargetingType> {

	typedef DensityMatrixBase<TargetingType> BaseType;
	typedef typename TargetingType::LeftRightSuperType LeftRightSuperType;
	typedef typename TargetingType::BasisWithOperatorsType BasisWithOperatorsType;
	typedef typename BasisWithOperatorsType::BasisType  BasisType;
	typedef typename BasisWithOperatorsType::SparseMatrixType SparseMatrixType;
	typedef typename TargetingType::TargetVectorType::value_type ComplexOrRealType;
	typedef typename BaseType::BlockDiagonalMatrixType BlockDiagonalMatrixType;
	typedef typename BasisType::FactorsType FactorsType;
	typedef typename PsimagLite::Real<ComplexOrRealType>::Type RealType;
	typedef typename DensityMatrixBase<TargetingType>::Params ParamsType;
	typedef typename BasisType::BlockType VectorSizeType;

public:

	typedef typename BlockDiagonalMatrixType::BuildingBlockType BuildingBlockType;

	DensityMatrixSu2(const TargetingType& target,
	                 const LeftRightSuperType& lrs,
	                 const ParamsType& p)
	    : pBasis_((p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.left() :
	                                                                              lrs.right()),
	      data_(pBasis_),
	      mMaximal_(data_.blocks()),
	      direction_(p.direction),
	      debug_(p.debug)
	{
		check();
		BuildingBlockType matrixBlock;

		const BasisWithOperatorsType& pBasisSummed =
		        (p.direction == ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) ? lrs.right() :
		                                                                        lrs.left();

		for (SizeType m = 0; m < pBasis_.partition() - 1; ++m) {
			// Definition: Given partition p with (j m)
			// findMaximalPartition(p) returns the partition p' (with j,j)

			if (BasisType::useSu2Symmetry()) {
				mMaximal_[m] = findMaximalPartition(m,pBasis_);
				//if (enforceSymmetry && SizeType(m)!=mMaximal_[m]) continue;
				// we'll fill non-maximal partitions later
			}

			SizeType bs = pBasis_.partition(m+1)-pBasis_.partition(m);

			matrixBlock.resize(bs, bs, static_cast<ComplexOrRealType>(0.0));

			for (SizeType i=pBasis_.partition(m);i<pBasis_.partition(m+1);i++) {
				for (SizeType j=pBasis_.partition(m);j<pBasis_.partition(m+1);j++) {

					matrixBlock(i-pBasis_.partition(m),j-pBasis_.partition(m))=
					        densityMatrixAux(i,j,target,pBasisSummed,lrs.super(),p.direction);

				}
			}

			data_.setBlock(m,pBasis_.partition(m),matrixBlock);
		}

		if (debug_) areAllMsEqual(pBasis_);
	}

	virtual const BlockDiagonalMatrixType& operator()()
	{
		return data_;
	}

	void diag(typename PsimagLite::Vector<RealType>::Type& eigs,char jobz)
	{
		DiagBlockDiagMatrix<BlockDiagonalMatrixType>::diagonalise(data_,eigs,jobz);

		//make sure non-maximals are equal to maximals
		// this is needed because otherwise there's no assure that m-independence
		// is achieved due to the non unique phase of eigenvectors of the density matrix
		for (SizeType m=0;m<data_.blocks();m++) {

			SizeType p = mMaximal_[m];
			if (SizeType(m)==p) continue; // we already did these ones

			data_.setBlock(m,data_.offsetsRows(m),data_(p));
		}

		if (debug_) areAllMsEqual(pBasis_);

		check2();
	}

	friend std::ostream& operator<<(std::ostream& os,
	                                const DensityMatrixSu2& dm)
	{
		for (SizeType m=0;m<dm.data_.blocks();m++) {
			// Definition: Given partition p with (j m)
			// findMaximalPartition(p) returns the partition p' (with j,j)
			std::pair<SizeType,SizeType> jm1 = dm.pBasis_.jmValue(dm.pBasis_.partition(m));
			SizeType ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<"partitionNumber="<<m<<" j="<<jm1.first;
			os<<" m= "<<jm1.second<<" ne="<<ne<<"\n";
			os<<dm.data_(m)<<"\n";
		}

		return os;
	}

private:

	void check()
	{
		if (!debug_) return;

		for (SizeType m = 0; m < data_.blocks(); ++m) {
			// Definition: Given partition p with (j m)
			// findMaximalPartition(p) returns the partition p' (with j,j)
			SizeType p = mMaximal_[m];
			if (m==p) continue;
			//is data_(m)==data_(p) ?
			check(m,data_(m),p,data_(p));
		}
	}

	void check2()
	{
		if (!debug_) return;
		for (SizeType m = 0; m < data_.blocks(); ++m) {
			// Definition: Given partition p with (j m)
			// findMaximalPartition(p) returns the partition p' (with j,j)
			SizeType p = mMaximal_[m];
			if (m==p) continue;
			//is data_(m)==data_(p) ?
			check2(m,data_(m),p,data_(p));

		}
	}

	SizeType findMaximalPartition(SizeType p, const BasisWithOperatorsType& pBasis)
	{
		VectorSizeType pBasisElectrons;
		pBasis.su2ElectronsBridge(pBasisElectrons);
		std::pair<SizeType,SizeType> jm2 = pBasis.jmValue(pBasis.partition(p));
		SizeType ne2 = pBasisElectrons[pBasis.partition(p)];
		if (jm2.first==jm2.second) return p;
		for (SizeType m=0;m<pBasis.partition()-1;m++) {
			std::pair<SizeType,SizeType> jm1 = pBasis.jmValue(pBasis.partition(m));
			SizeType ne1 = pBasisElectrons[pBasis.partition(m)];
			if (jm1.first==jm2.first && jm1.first==jm1.second && ne1==ne2) return m;
		}

		throw PsimagLite::RuntimeError("findMaximalPartition : none found\n");
	}

	//! only used for debugging
	bool areAllMsEqual(const BasisWithOperatorsType&)
	{
		return true;
	}

	ComplexOrRealType densityMatrixAux(SizeType alpha1,
	                                          SizeType alpha2,
	                                          const TargetingType& target,
	                                          const BasisWithOperatorsType& pBasisSummed,
	                                          const BasisType& pSE,
	                                          ProgramGlobals::DirectionEnum direction)
	{
		ComplexOrRealType sum=0;
		// The g.s. has to be treated separately because it's
		// usually a vector of RealType, whereas
		// the other targets might be complex,
		// and C++ generic programming capabilities are weak... we need D!!!
		if (target.includeGroundStage())
			sum +=  densityMatrixHasFactors(alpha1,
			                                alpha2,
			                                target.gs(),
			                                pBasisSummed,
			                                pSE,
			                                direction)*target.gsWeight();

		for (SizeType i = 0; i < target.size(); ++i)
			sum += densityMatrixHasFactors(alpha1,
			                               alpha2,
			                               target(i),
			                               pBasisSummed,
			                               pSE,
			                               direction)*target.weight(i)/target.normSquared(i);

		return sum;
	}

	template<typename TargetVectorType>
	ComplexOrRealType densityMatrixHasFactors(SizeType alpha1,
	                                                 SizeType alpha2,
	                                                 const TargetVectorType& v,
	                                                 const BasisWithOperatorsType& pBasisSummed,
	                                                 const BasisType& pSE,
	                                                 ProgramGlobals::DirectionEnum direction)
	{
		int ne = pBasisSummed.size();
		int ns = pSE.size()/ne;
		SizeType total=pBasisSummed.size();
		if (direction != ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
			ns=pBasisSummed.size();
			ne=pSE.size()/ns;
		}

		ComplexOrRealType sum=0;

		// Make sure we don't copy just get the reference here!!
		const FactorsType* fptr = pSE.getFactors();
		assert(fptr);
		const FactorsType& factors = *fptr;

		for (SizeType beta=0;beta<total;beta++) {
			// sum over environ:
			int i1 = alpha1+beta*ns;
			int i2 = alpha2+beta*ns;
			// sum over system:
			if (direction != ProgramGlobals::DirectionEnum::EXPAND_SYSTEM) {
				i1 = beta + alpha1*ns;
				i2 = beta + alpha2*ns;
			}
			for (int k1=factors.getRowPtr(i1);k1<factors.getRowPtr(i1+1);k1++) {
				int eta1 = factors.getCol(k1);
				int ii = pSE.permutationInverse(eta1);
				for (int k2=factors.getRowPtr(i2);k2<factors.getRowPtr(i2+1);k2++) {
					int eta2 =  factors.getCol(k2);
					int jj = pSE.permutationInverse(eta2);

					ComplexOrRealType tmp3= v.slowAccess(ii)*
					        PsimagLite::conj(v.slowAccess(jj)) *
					        factors.getValue(k1) * factors.getValue(k2);
					sum += tmp3;
				}
			}
		}
		return sum;
	}

	//! only used for debugging
	void check(SizeType p1,
	           const BuildingBlockType& bp1,
	           SizeType p2,
	           const BuildingBlockType& bp2)
	{
		if (bp1.rows()!=bp2.rows()) {
			std::cerr<<"row size different "<<bp1.rows()<<"!="<<bp2.rows()<<"\n";
			std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
			throw PsimagLite::RuntimeError("Density Matrix Check: failed\n");
		}

		if (bp1.cols()!=bp2.cols()) {
			std::cerr<<"col size different "<<bp1.cols()<<"!="<<bp2.cols()<<"\n";
			std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
			throw PsimagLite::RuntimeError("Density Matrix Check: failed\n");
		}

		if (!debug_) return;

		for (SizeType i=0;i<bp1.rows();i++) {
			for (SizeType j=0;j<bp1.cols();j++) {
				RealType x = PsimagLite::norm(bp1(i,j)-bp2(i,j));
				if (x>1e-5) {
					std::cerr<<bp1(i,j)<<"!="<<bp2(i,j)<<" i="<<i<<" j= "<<j<<"\n";
					std::cerr<<"difference="<<x<<" p1="<<p1<<" p2="<<p2<<"\n";
					std::cerr<<data_(p1)<<"\n";
					std::cerr<<"******************\n";
					std::cerr<<data_(p2)<<"\n";
					throw PsimagLite::RuntimeError("Density Matrix Check: failed (differ)\n");
				}
			}
		}
	}

	//! only used for debugging
	void check2(SizeType p1,
	            const BuildingBlockType& bp1,
	            SizeType p2,
	            const BuildingBlockType& bp2)
	{
		ComplexOrRealType alpha=1.0,beta=0.0;
		int n =bp1.cols();
		PsimagLite::Matrix<ComplexOrRealType> result(n,n);
		psimag::BLAS::GEMM('C',
		                   'N',
		                   n,
		                   n,
		                   n,
		                   alpha,
		                   &(bp1(0,0)),
		                   n,
		                   &(bp2(0,0)),
		                   n,
		                   beta,
		                   &(result(0,0)),
		                   n);
		if (!PsimagLite::isTheIdentity(result)) {
			//utils::matrixPrint(result,std::cerr);
			std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
			throw PsimagLite::RuntimeError("Density Matrix Check2: failed\n");
		}
	}

	const BasisWithOperatorsType& pBasis_;
	BlockDiagonalMatrixType data_;
	VectorSizeType mMaximal_;
	ProgramGlobals::DirectionEnum direction_;
	bool debug_;
}; // class DensityMatrixSu2
} // namespace Dmrg

/*@}*/
#endif
