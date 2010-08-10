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

/*! \file DensityMatrixSu2.h
 *
 *  
 *
 */
#ifndef DENSITY_MATRIX_SU2_H
#define DENSITY_MATRIX_SU2_H

#include "Utils.h"
#include "BlockMatrix.h"
#include "DensityMatrixBase.h"

namespace Dmrg {
	//!
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
		typename TargettingType
		>
	class DensityMatrixSu2 : public DensityMatrixBase<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType> {
		typedef typename DmrgBasisWithOperatorsType::SparseMatrixType SparseMatrixType;
		typedef typename TargettingType::TargetVectorType::value_type DensityMatrixElementType;
		typedef BlockMatrix<DensityMatrixElementType,psimag::Matrix<DensityMatrixElementType> > BlockMatrixType;
		typedef typename DmrgBasisType::FactorsType FactorsType;
		enum {EXPAND_SYSTEM = TargettingType::EXPAND_SYSTEM };
		
	public:
		typedef typename BlockMatrixType::BuildingBlockType BuildingBlockType;
		
		DensityMatrixSu2(
			const TargettingType& target,
			const DmrgBasisWithOperatorsType& pBasis,
			const DmrgBasisWithOperatorsType& pBasisSummed,
			const DmrgBasisType& pSE,
			size_t direction,bool debug=false,bool verbose=false) : data_(pBasis.size() ,pBasis.partition()-1),mMaximal_(pBasis.partition()-1),pBasis_(pBasis),
				debug_(debug),verbose_(verbose)
		{
		}
		
		BlockMatrixType& operator()()
		{
			return data_;
		}
		
		size_t rank() { return data_.rank(); }
		
		void check(int direction)
		{
			if (!debug_) return;
			
			if (verbose_) std::cerr<<"CHECKING DENSITY-MATRIX WITH OPTION="<<direction<<"\n";
			for (size_t m=0;m<data_.blocks();m++) {
				// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
				size_t p = mMaximal_[m];
				if (m==p) continue;
				//is data_(m)==data_(p) ?
				check(m,data_(m),p,data_(p));
				
			}
		}

		void check2(int direction)
		{
			if (!debug_) return;
			if (verbose_) std::cerr<<"CHECKING DMRG-TRANFORM WITH OPTION="<<direction<<"\n";
			for (size_t m=0;m<data_.blocks();m++) {
				// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
				size_t p = mMaximal_[m];
				if (m==p) continue;
				//is data_(m)==data_(p) ?
				check2(m,data_(m),p,data_(p));
				
			}
		}

		template<typename ConcurrencyType>
		void diag(std::vector<RealType>& eigs,char jobz,ConcurrencyType& concurrency)
		{
			diagonalise<DensityMatrixElementType,RealType,ConcurrencyType>(data_,eigs,jobz,concurrency);
			
			//make sure non-maximals are equal to maximals
			// this is needed because otherwise there's no assure that m-independence
			// is achieved due to the non unique phase of eigenvectors of the density matrix
			for (size_t m=0;m<data_.blocks();m++) {
				
				size_t p = mMaximal_[m];
				if (size_t(m)==p) continue; // we already did these ones

				data_.setBlock(m,data_.offsets(m),data_(p));
			}
			if (verbose_) std::cerr<<"After diagonalise\n";

			if (debug_) areAllMsEqual(pBasis_);
			
		}
		
		void init(
				const TargettingType& target,
				DmrgBasisWithOperatorsType const &pBasis,
				const DmrgBasisWithOperatorsType& pBasisSummed,
				DmrgBasisType const &pSE,
				int direction)
		{
			BuildingBlockType matrixBlock;
			
			for (size_t m=0;m<pBasis.partition()-1;m++) {
				// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
				
				if (DmrgBasisType::useSu2Symmetry()) {
					mMaximal_[m] = findMaximalPartition(m,pBasis);
					//if (enforceSymmetry && size_t(m)!=mMaximal_[m]) continue; 
					// we'll fill non-maximal partitions later
				}

				size_t bs = pBasis.partition(m+1)-pBasis.partition(m);
				
				matrixBlock.resize(bs,bs);
				
				for (size_t i=pBasis.partition(m);i<pBasis.partition(m+1);i++) {
					for (size_t j=pBasis.partition(m);j<pBasis.partition(m+1);j++) {
						
						matrixBlock(i-pBasis.partition(m),j-pBasis.partition(m))=
							densityMatrixAux(i,j,target,pBasisSummed,pSE,direction);
						
					}
				}
				data_.setBlock(m,pBasis.partition(m),matrixBlock);
			}

			if (verbose_) {
				std::cerr<<"DENSITYMATRIXPRINT option="<<direction<<"\n";
				std::cerr<<(*this);
				std::cerr<<"***********\n";
				std::cerr<<"Calling ae from init()...\n";
			}
			if (debug_) areAllMsEqual(pBasis);
		}

		template<
			typename RealType_,
			typename DmrgBasisType_,
			typename DmrgBasisWithOperatorsType_,
   			typename TargettingType_
			> 
		friend std::ostream& operator<<(std::ostream& os,
				const DensityMatrixSu2<RealType_,
    					DmrgBasisType_,DmrgBasisWithOperatorsType_,TargettingType_>& dm);
	private:
		BlockMatrixType data_;
		std::vector<size_t> mMaximal_;
		const DmrgBasisWithOperatorsType& pBasis_;
		bool debug_,verbose_;
		
		size_t findMaximalPartition(size_t p,DmrgBasisWithOperatorsType const &pBasis)
		{
			std::pair<size_t,size_t> jm2 = pBasis.jmValue(pBasis.partition(p));
			size_t ne2 = pBasis.electrons(pBasis.partition(p));
			if (jm2.first==jm2.second) return p;
			for (size_t m=0;m<pBasis.partition()-1;m++) { 
				std::pair<size_t,size_t> jm1 = pBasis.jmValue(pBasis.partition(m));
				size_t ne1 = pBasis.electrons(pBasis.partition(m));
				if (jm1.first==jm2.first && jm1.first==jm1.second && ne1==ne2) return m;
			}
			throw std::runtime_error("	findMaximalPartition : none found\n");
		}

		//! only used for debugging
		bool areAllMsEqual(DmrgBasisWithOperatorsType const &pBasis)
		{
			for (size_t m=0;m<pBasis.partition()-1;m++) {
				std::pair<size_t,size_t> jmPair1 =  pBasis.jmValue(pBasis.partition(m));
				size_t ne1 = pBasis.electrons(pBasis.partition(m));
				for (size_t p=0;p<pBasis.partition()-1;p++) {
					std::pair<size_t,size_t> jmPair2 =  pBasis.jmValue(pBasis.partition(p));
					size_t ne2 = pBasis.electrons(pBasis.partition(p));

					if (jmPair1.first == jmPair2.first && ne1==ne2) {
						
						if (!almostEqual(data_(m),data_(p),1e-5)) {
							std::cerr<<"Checking m="<<m<<" against p="<<p<<"\n";
							throw std::runtime_error("error\n");
						}
					}
				}
			}
			return true;
		}
		DensityMatrixElementType densityMatrixAux(size_t alpha1,size_t alpha2,const TargettingType& target,
			DmrgBasisWithOperatorsType const &pBasisSummed,DmrgBasisType const &pSE,size_t direction)
		{
			DensityMatrixElementType sum=0;
			// The g.s. has to be treated separately because it's usually a vector of double, whereas
			// the other targets might be complex, and C++ generic programming capabilities are weak... we need D!!!
			if (target.includeGroundStage()) sum +=  densityMatrixHasFactors(alpha1,alpha2,target.gs(),
			    		pBasisSummed,pSE,direction)*target.gsWeight();
			
			for (size_t i=0;i<target.size();i++) 
				sum += densityMatrixHasFactors(alpha1,alpha2,target(i),
					pBasisSummed,pSE,direction)*target.weight(i)/target.normSquared(i);

			return sum;
		}
		
		template<typename TargetVectorType>
		DensityMatrixElementType densityMatrixHasFactors(size_t alpha1,size_t alpha2,const TargetVectorType& v,
			DmrgBasisWithOperatorsType const &pBasisSummed,DmrgBasisType const &pSE,size_t direction)
		{
			int ne = pBasisSummed.size();
			int ns = pSE.size()/ne;
			size_t total=pBasisSummed.size();
			if (direction!=EXPAND_SYSTEM) {
				ns=pBasisSummed.size();
				ne=pSE.size()/ns;
			}
				
			DensityMatrixElementType sum=0;

			// Make sure we don't copy just get the reference here!!
			const FactorsType& factors = pSE.getFactors();

			for (size_t beta=0;beta<total;beta++) {
				// sum over environ:
				int i1 = alpha1+beta*ns;
				int i2 = alpha2+beta*ns;
				// sum over system:
				if (direction!=EXPAND_SYSTEM) {
					i1 = beta + alpha1*ns;
					i2 = beta + alpha2*ns;
				}		
				for (int k1=factors.getRowPtr(i1);k1<factors.getRowPtr(i1+1);k1++) {
					int eta1 = factors.getCol(k1);
					int ii = pSE.permutationInverse(eta1);
					for (int k2=factors.getRowPtr(i2);k2<factors.getRowPtr(i2+1);k2++) {
						int eta2 =  factors.getCol(k2);
						int jj = pSE.permutationInverse(eta2);
					
						DensityMatrixElementType tmp3= v[ii] * std::conj(v[jj]) * 
							factors.getValue(k1) * factors.getValue(k2);
						sum += tmp3;
					}
				}
			}
			return sum;
		}

		//! only used for debugging
		void check(size_t p1,const BuildingBlockType& bp1,size_t p2,const BuildingBlockType& bp2)
		{
			if (bp1.n_row()!=bp2.n_row()) {
				std::cerr<<"row size different "<<bp1.n_row()<<"!="<<bp2.n_row()<<"\n";
				std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
				std::cerr<<data_<<"\n";
				throw std::runtime_error("Density Matrix Check: failed\n");
			}
			if (bp1.n_col()!=bp2.n_col()) {
				std::cerr<<"col size different "<<bp1.n_col()<<"!="<<bp2.n_col()<<"\n";
				std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
				std::cerr<<data_<<"\n";
				throw std::runtime_error("Density Matrix Check: failed\n");
			}
			if (!debug_) return;

			for (size_t i=0;i<bp1.n_row();i++) {
				for (size_t j=0;j<bp1.n_col();j++) {
					RealType x =psimag::norm(bp1(i,j)-bp2(i,j)); 
					if (x>1e-5) {
						std::cerr<<bp1(i,j)<<"!="<<bp2(i,j)<<" i="<<i<<" j= "<<j<<"\n";
						std::cerr<<"difference="<<x<<" p1="<<p1<<" p2="<<p2<<"\n";
						std::cerr<<data_(p1)<<"\n";
						std::cerr<<"******************\n";
						std::cerr<<data_(p2)<<"\n";
						throw std::runtime_error("Density Matrix Check: failed (differ)\n");
					}
				}
			}
		}
		
		//! only used for debugging
		void check2(size_t p1,const BuildingBlockType& bp1,size_t p2,const BuildingBlockType& bp2)
		{ 
			DensityMatrixElementType alpha=1.0,beta=0.0;
			int n =bp1.n_col();
			psimag::Matrix<DensityMatrixElementType> result(n,n);
			psimag::BLAS::GEMM('C','N',n,n,n,alpha,&(bp1(0,0)),n,&(bp2(0,0)),n,beta,&(result(0,0)),n);
			if (!psimag::isTheIdentity<DensityMatrixElementType,RealType>(result)) {
				//utils::matrixPrint(result,std::cerr);
				std::cerr<<"p1="<<p1<<" p2="<<p2<<"\n";
				throw std::runtime_error("Density Matrix Check2: failed\n");
			}
					
		}
		
	}; // class DensityMatrixSu2
	
	template<
		typename RealType,
		typename DmrgBasisType,
		typename DmrgBasisWithOperatorsType,
  		typename TargettingType
		> 
	std::ostream& operator<<(std::ostream& os,
				const DensityMatrixSu2<RealType,DmrgBasisType,DmrgBasisWithOperatorsType,TargettingType>& dm)
	{
		//std::cerr<<"PRINTING DENSITY-MATRIX WITH OPTION="<<option<<"\n";
		for (size_t m=0;m<dm.data_.blocks();m++) {
			// Definition: Given partition p with (j m) findMaximalPartition(p) returns the partition p' (with j,j)
			std::pair<size_t,size_t> jm1 = dm.pBasis_.jmValue(dm.pBasis_.partition(m));
			size_t ne = dm.pBasis_.electrons(dm.pBasis_.partition(m));
			os<<"partitionNumber="<<m<<" j="<<jm1.first<<" m= "<<jm1.second<<" ne="<<ne<<"\n"; 
			os<<dm.data_(m)<<"\n";
		}
		return os;
	}
} // namespace Dmrg

/*@}*/
#endif

 
